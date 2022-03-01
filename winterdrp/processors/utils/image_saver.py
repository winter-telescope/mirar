import os
from winterdrp.paths import output_path, output_dir
from winterdrp.processors.base_processor import BaseProcessor

latest_save_key = "SAVEPATH"


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

        try:
            os.makedirs(output_dir(dir_root=self.save_dir, sub_dir=sub_dir))
        except OSError:
            pass

        for i, img in enumerate(images):

            header = headers[i]

            path = output_path(header["BASENAME"], dir_root=self.save_dir, sub_dir=sub_dir)

            self.save_fits(img, header, path)
            header[latest_save_key] = path

        return images, headers
