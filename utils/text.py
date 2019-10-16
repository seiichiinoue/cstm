import argparse
import os
import MeCab

class Text(object):
    def __init__(self, path="data/wiki.txt"):
        self.wakati = MeCab.Tagger('-Owakati')
        self.chasen = MeCab.Tagger('-Ochasen')
        with open(path, "r", encoding="utf-8") as f:
            self.text = f.readlines()
        for i in range(len(self.text)):
            self.text[i] = self.text[i].replace("\n", "")

    def _wakati_ja(self, path):
        wakatied = []
        for i in range(len(self.text)):
            new_t = self.wakati.parse(self.text[i]).strip("\n")
            wakatied.append(new_t)
        
        with open(path, "w") as f:
            for t in wakatied:
                f.write(t+'\n')
        return None
    
    def _wakati_en(self, path):
        wakatied = []
        for i in range(len(self.text)):
            new_t = self.text[i].split(" ").strip("\n")
            wakatied.append(new_t)
        with open(path, "w") as f:
            for t in wakatied:
                f.write(t+'\n')
        return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='this script for text processing.')
    parser.add_argument('--lang_type', help='english or japanese')
    parser.add_argument('--tar_path', help='directory located raw file')
    parser.add_argument('--save_path', help='directory save processed text')
    args = parser.parse_args()
    filelist = os.listdir(args.tar_path)
    if args.lang_type.lower() == 'japanese':
        for file in filelist:
            if file.endswith(".txt"):
                t = Text(args.tar_path+file)
                t._wakati_ja(args.save_path+file)
    elif args.lang_type.lower() == 'english':
        for file in filelist:
            if file.endswith(".txt"):
                t = Text(args.tar_path+file)
                t._wakati_en(args.save_path+file)
    else:
        print("Error: this script does not support language of got file")