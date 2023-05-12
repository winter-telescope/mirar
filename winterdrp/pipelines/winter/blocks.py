from winterdrp.pipelines.winter.generator import winter_reference_generator
from winterdrp.processors.reference import GetReferenceImage
from winterdrp.processors.utils import ImageDebatcher, ImageSaver

refbuild = [
    ImageDebatcher(),
    GetReferenceImage(
        ref_image_generator=winter_reference_generator,
    ),
    ImageSaver(output_dir_name="stacked_ref"),
]
