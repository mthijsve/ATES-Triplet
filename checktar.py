import tarfile
import os
import shutil

# Define the path to the tar.gz file and the destination directory
tar_file = 'output/10e12b.tar.gz'
dest_dir = 'Paper1 output'

# Create the destination directory if it doesn't exist
if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)

# Open the tar.gz archive for reading
with tarfile.open(tar_file, 'r:gz') as archive:
    # Loop through all files in the archive
    for member in archive.getmembers():
        # Check if the file has a .csv extension
        if member.name.endswith('.csv'):
            # Extract the file to the destination directory
            archive.extract(member, dest_dir)
print('done', tar_file)