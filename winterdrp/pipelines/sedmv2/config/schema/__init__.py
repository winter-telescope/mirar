import os

sedmv2_schema_dir = os.path.dirname(__file__)


def get_sedmv2_schema_path(
        db_name: str
) -> str:
    return os.path.join(sedmv2_schema_dir, f"{db_name}.sql")
