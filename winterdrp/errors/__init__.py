from winterdrp.paths import base_name_key
from astropy.time import Time


class ProcessorError(BaseException):
    pass


class ErrorReport:

    def __init__(
            self,
            error,
            processor_name,
            contents: list[str]
    ):
        self.error = error
        self.processor_name = processor_name
        self.contents = contents
        self.t_error = Time.now()

    def generate_log_message(self):
        return f"Error for processor {self.processor_name} at time {self.t_error} UT: " \
               f"{type(self.error).__name__} affected batch of length {len(self.contents)}."