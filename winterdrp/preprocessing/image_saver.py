from winterdrp.preprocessing.base_processor import BaseProcessor


class ImageSaver(BaseProcessor):

    base_key = "save"

    def _apply_to_images(
            self,
            images: list,
            headers: list,
            sub_dir: str = ""
    ) -> (list, list):
        print("run")
        return images, headers
