import tarfile

tar = tarfile.open("output.tar.gz", "r:gz")
tar.getmembers()