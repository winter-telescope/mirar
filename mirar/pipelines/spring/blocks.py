# TODO : write blocks for spring pipeline
from mirar.pipelines.spring.load_spring_image import load_raw_spring_image
from mirar.processors.utils.image_loader import ImageLoader


load_raw = [
    ImageLoader(
    input_sub_dir="raw",
    load_image=load_raw_spring_image)
]

csvlog = [

    ## TODO : Add logging block here ##
]