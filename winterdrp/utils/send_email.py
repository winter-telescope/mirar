import os
import smtplib
from email.message import EmailMessage
import ssl
import getpass

port = 465  # For SSL


def send_email(
        email_recipients: str | list[str],
        email_subject: str,
        email_text: str,
        email_sender: str,
        email_password: str = os.getenv("EMAIL_PASSWORD")
):

    # Create a text/plain message
    msg = EmailMessage()
    msg.set_content(email_text)

    if not isinstance(email_recipients, list):
        email_recipients = [email_recipients]

    # # me == the sender's email address
    # # you == the recipient's email address
    msg['Subject'] = email_subject
    msg['From'] = email_sender
    msg['To'] = ', '.join(email_recipients)

    print(msg)

    if email_password is None:
        email_password = getpass.getpass()

    # Create a secure SSL context
    context = ssl.create_default_context()

    with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
        server.login(email_sender, email_password)
        server.send_message(msg)


if __name__ == "__main__":

    send_email(
        email_sender="winter.data.reduction.pipeline@gmail.com",
        email_recipients="rdstein@caltech.edu",
        email_subject="Tester",
        email_text="a"
    )
