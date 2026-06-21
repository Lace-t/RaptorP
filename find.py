import os
import sys
import re

def main():
    keyword = str(sys.argv[1])
    target_dir = sys.argv[2]
    pattern = re.compile(r'\b' + keyword + r'\b')
    
    # Optional: Output results to a file
    #output_file = open('files_with_'+keyword+'.txt', 'w')
    
    for root, dirs, files in os.walk(target_dir):
        for file in files:
            if file.endswith(("c",'.h')):
                file_path = os.path.join(root, file)
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                        if pattern.search(content):
                            print(f"Found '{keyword}' in: {file_path}")
                            #output_file.write(file_path + '\n')
                except UnicodeDecodeError:
                    print(f"Skipping file due to decoding error: {file_path}")
                except Exception as e:
                    print(f"Error reading file {file_path}: {str(e)}")
    
    #output_file.close()

if __name__ == "__main__":
    main()
