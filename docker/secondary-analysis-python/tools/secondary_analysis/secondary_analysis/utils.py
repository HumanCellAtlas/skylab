import requests
import logging
import time


def get_file_by_uuid(file_id, dss_url):
    """
    Retrieve a file from the Human Cell Atlas data storage service by its id.
    :param str file_id: The id of the file to retrieve
    :param str dss_url: The url for the HCA data storage service, e.g. "https://dss.staging.data.humancellatlas.org/v1"
    :return: dict file contents
    """
    url = '{dss_url}/files/{file_id}?replica=gcp'.format(
        dss_url=dss_url, file_id=file_id)
    logging.info('GET {0}'.format(url))
    response = requests.get(url)
    logging.info(response.status_code)
    logging.info(response.text)
    return response.json()


def get_manifest_files(bundle_uuid, bundle_version, dss_url, timeout_seconds, retry_seconds):
    """
    Retrieve manifest.json file for a given bundle uuid and version.
    :param str bundle_uuid: Bundle unique id
    :param str bundle_version: Timestamp of bundle creation, e.g. "2017-10-23T17:50:26.894Z"
    :param str dss_url: The url for the Human Cell Atlas data storage service, e.g. "https://dss.staging.data.humancellatlas.org/v1"
    :param int timeout_seconds: Seconds before allowing the request to timeout
    :param int retry_seconds: Seconds between retrying the request to get the manifest file
    :return: {
                'name_to_meta': dict mapping <str file name>: <dict file metadata>,
                'url_to_name': dict mapping <str file url>: <str file name>
             }
    """
    url = '{dss_url}/bundles/{bundle_uuid}?version={bundle_version}&replica=gcp&directurls=true'.format(
        dss_url=dss_url, bundle_uuid=bundle_uuid, bundle_version=bundle_version)
    start = time.time()
    current = start
    # Retry in a loop because of intermittent 5xx errors from dss
    while current - start < timeout_seconds:
        logging.info('GET {0}'.format(url))
        response = requests.get(url)
        logging.info(response.status_code)
        logging.info(response.text)
        if 200 <= response.status_code <= 299:
            break
        time.sleep(retry_seconds)
        current = time.time()
    manifest = response.json()

    bundle = manifest['bundle']
    name_to_meta = {}
    url_to_name = {}
    for f in bundle['files']:
        name_to_meta[f['name']] = f
        url_to_name[f['url']] = f['name']

    return {
        'name_to_meta': name_to_meta,
        'url_to_name': url_to_name
    }
