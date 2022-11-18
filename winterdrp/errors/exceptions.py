class BaseProcessorError(Exception):
    pass


class ProcessorError(BaseProcessorError):
    pass


class NoncriticalProcessingError(BaseProcessorError):
    pass


class ImageNotFoundError(ProcessorError, FileNotFoundError):
    pass
