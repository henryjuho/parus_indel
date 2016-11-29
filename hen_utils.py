def file_location(file_string):
    """
    function that strip the filename off a path

    :param file_string: str
    :return: path str
    """
    return file_string[:file_string.rfind('/')+1]


def extract_filename(file_string):
    """
    function that strips a path of a filename

    :param file_string: str
    :return: file str
    """
    return file_string[file_string.rfind('/')+1:]