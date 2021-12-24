from winterdrp.paths import output_dir
from winterdrp.preprocessing.base_processor import BaseProcessor


class ImageSaver(BaseProcessor):

    base_key = "save"

    def __init__(
            self,
            instrument_vars: dict,
            save_dir: str,
            *args,
            **kwargs
    ):
        super().__init__(instrument_vars, *args, **kwargs)
        self.save_dir = save_dir

    def _apply_to_images(
            self,
            images: list,
            headers: list,
            sub_dir: str = ""
    ) -> (list, list):

        for i, img in enumerate(images):

            paths = output_dir(self.save_dir, sub_dir=sub_dir)

            self.save_fits(img, headers[i], path)

        return images, headers
