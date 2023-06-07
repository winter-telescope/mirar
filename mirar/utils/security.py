"""
Util functions for password generation
"""
import secrets
import string


def generate_key(length: int = 20) -> str:
    """
    Generate an alphanumeric password of length N, from
    https://docs.python.org/3/library/secrets.html#recipes-and-best-practices

    :param length: length
    :return: Password
    """
    alphabet = string.ascii_letters + string.digits
    return "".join(secrets.choice(alphabet) for _ in range(length))
