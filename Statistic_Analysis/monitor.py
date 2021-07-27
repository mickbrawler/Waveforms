import os
import time
import argparse

def monitor(directory, dir_length):

    while True: 
        results = os.listdir(directory)
        time.sleep(1)
        if len(results) == (dir_length):
            break

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type=str)
    parser.add_argument('--dir_length', type=int)
    args = parser.parse_args()

    monitor(args.directory,args.dir_length)
