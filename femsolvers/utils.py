import os


def fix_docker_permissions(root_dir: str):
    """
    Fix docker UID issue. NOTE: One should use this very carefully.
    This function gives permissions for EVERYONE to do EVERYTHING with all subdirs of
    root_dir.

    :param str root_dir:
    :return:
    """
    for root, dirs, files in os.walk(root_dir):
        do_what_you_want = 0o777

        os.chmod(root, do_what_you_want)

        for directory in dirs:
            os.chmod(os.path.join(root, directory), do_what_you_want)

        for file in files:
            os.chmod(os.path.join(root, file), do_what_you_want)