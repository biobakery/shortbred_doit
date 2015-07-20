#! usr/bin/env python

"""
Map the Broad folder IDs to the more common
HPFS donor IDs and sample metadata, e.g. DNA vs. RNA, body site, etc.
"""

import glob
import os

# headers
print "\t".join(["headers", "folder", "donor", "method", "site", "type"])

# process folders
root = "hpfs/broad/B452"
folders = glob.glob(os.path.join(root, "*"))
for folder in folders:
    if os.path.isdir(folder):
        foldername = os.path.split(folder)[-1]
        subfolders = glob.glob(os.path.join(folder, "*"))
        for subfolder in subfolders:
            if os.path.isdir(subfolder):
                fh = open(os.path.join(subfolder, "metadata.txt"))
                lines = fh.readlines()
                fh.close()
                items = lines[4].strip().split()
                sampleid = items[-1]
                (donor, method) = sampleid.split("_")
                nttype = "DNA" if foldername[0] == "G" else "RNA"
                site = "Saliva" if method == "Saliva" else "Stool"
                method = "Whole" if method == "Saliva" else method
                print "\t".join([foldername, subfolder, donor, method, site, nttype]) 
