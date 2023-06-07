"""
Module containing gmail integration functions
"""
import getpass
import gzip
import logging
import os
import smtplib
import ssl
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

GMAIL_PORT = 465  # For SSL


def send_gmail(
    email_recipients: str | list[str],
    email_subject: str,
    email_text: str,
    email_sender: str = os.getenv("WATCHDOG_EMAIL"),
    email_password: str = os.getenv("WATCHDOG_EMAIL_PASSWORD"),
    attachments: Optional[str | list[str]] = None,
    auto_compress: bool = True,
):
    """
    Function to send an email to a list of recipients from a gmail account.

    :param email_recipients: recipients for email
    :param email_subject: subject for the email
    :param email_text: Text to send
    :param email_sender: Gmail to send from
    :param email_password: Password for sender gmail account
    :param attachments: Any files to attach
    :param auto_compress: Boolean to compress large attachments before sending
    :return:
    """
    # pylint: disable=too-many-arguments

    # Create a text/plain message
    msg = MIMEMultipart()
    msg.attach(MIMEText(email_text))

    if not isinstance(email_recipients, list):
        email_recipients = [email_recipients]

    msg["Subject"] = email_subject
    msg["From"] = email_sender
    msg["To"] = ", ".join(email_recipients)

    if attachments is None:
        attachments = []

    if not isinstance(attachments, list):
        attachments = [attachments]

    for file_path in attachments:
        if os.path.exists(file_path):
            base_name = os.path.basename(file_path)

            if not isinstance(file_path, Path):
                file_path = Path(file_path)

            with open(file_path, "rb") as attachment:
                if np.logical_and(
                    auto_compress, file_path.stat().st_size > (1024 * 1024)
                ):
                    data = gzip.compress(attachment.read())
                    base_name += ".gzip"
                else:
                    data = attachment.read()

                part = MIMEApplication(data, Name=base_name)

            # After the file is closed
            part["Content-Disposition"] = f"attachment; filename={base_name}"
            msg.attach(part)

        else:
            logger.warning(f"Attachment {file_path} not found, skipping.")

    logger.info(f"Sending email to {email_recipients}")

    if email_password is None:
        email_password = getpass.getpass()

    # Create a secure SSL context
    context = ssl.create_default_context()

    with smtplib.SMTP_SSL("smtp.gmail.com", GMAIL_PORT, context=context) as server:
        server.login(email_sender, email_password)
        server.send_message(msg)
