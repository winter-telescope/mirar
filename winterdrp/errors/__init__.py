from winterdrp.paths import base_name_key

class ErrorReport:

    def __init__(
            self,
            error,
            processor_name,
            images,
            headers
    ):
        self.error = error
        self.processor_name = processor_name
        self.filenames = [h[base_name_key] for h in headers]

    def generate_log_message(self):
        return "You messed up somewhere :)"