## foto-import.py
# -*- coding: utf-8 -*-

import os
import shutil
from datetime import datetime
from dateutil.relativedelta import relativedelta

#--------------------------------------------------------
# Input
src_dir = r'D:\Fotos und Filme\iCloud Photos'  
dst_dir = r'D:\Fotos und Filme\Fotos' 
delete = False
#--------------------------------------------------------

def copy_and_delete_files_by_date(src_dir, dst_dir):
    """
    Copy files from src_base_dir to dst_base_dir based on the date, then delete the original files.

    :param src_base_dir: Source directory where files are located.
    :param dst_base_dir: Destination base directory to copy files to.
    """
    if not os.path.exists(src_dir):
        print(f"Source directory '{src_dir}' does not exist.")
        return

    if not os.path.exists(dst_dir):
        print(f"Destination base directory '{dst_dir}' does not exist.")
        return

    count = 0
    # List files in the source directory
    for filename in os.listdir(src_dir):
        file_path = os.path.join(src_dir, filename)
        print(f"Checking file: {file_path}")

        if os.path.isfile(file_path):
            # Get the file creation and modification times
            try:
                file_date_ = os.path.getmtime(file_path)
                
                file_date = datetime.fromtimestamp(file_date_)
                print(f"File last modified date: {file_date}")
            except Exception as e:
                print(f"Error processing file {file_path}: {e}")
        else:
            print(f"{file_path} is not a file.")

        target_year_dir = os.path.join(dst_dir,file_date.strftime("%Y"))
        target_month_dir = os.path.join(target_year_dir, file_date.strftime("%m-%Y"))
        if not os.path.exists(target_year_dir):
            print(f"Creating target directory: {target_year_dir}")
            os.makedirs(target_year_dir)
        if not os.path.exists(target_month_dir):
            print(f"Creating target directory: {target_month_dir}")
            os.makedirs(target_month_dir)

        # Copy the file
        print(f"Copying {filename} to {target_month_dir}")
        shutil.copy2(file_path, target_month_dir)
            # Delete the original file
        if delete == True:
            os.remove(file_path)
            print(f"Copied and deleted {filename}")
        count += 1
        #print("the current iterable is: ", filename)
    print(f"Finished copying {count} file/s to its/their target directory/ies.")

copy_and_delete_files_by_date(src_dir, dst_dir)