# todo finish this later

task DownloadFiles {
  Array[File] links
  String local_storage_directory

  command {
    for link in ${links}; do
      if [[ $link == "gs://"* ]]; then
        # download
      elif [[ $link == "s3://"* ]]; then
        # download
      elif [[ $link == "https://"* ]] || [[ $link == "http://"* ]] || [[ $link == "ftp://"* ]]; then
        # download
      else
        # report errors
    done
  }

}