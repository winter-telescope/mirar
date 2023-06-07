from mirar.pipelines.winter.generator import winter_reference_generator
from mirar.processors.reference import GetReferenceImage
from mirar.processors.utils import ImageDebatcher, ImageSaver

refbuild = [
    ImageDebatcher(),
    GetReferenceImage(
        ref_image_generator=winter_reference_generator,
    ),
    ImageSaver(output_dir_name="stacked_ref"),
]
