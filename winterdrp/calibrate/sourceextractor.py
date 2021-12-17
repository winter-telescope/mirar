from winterdrp.paths import sextractor_path, get_docker_client, docker_image

client = get_docker_client()

print(client.images.list())