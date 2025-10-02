'''
This script converts gff3 to BED format for gene features

accepts two arguements:
    first argument is the input file
    second argument is the output file

the script will parse each line, separating each column into local directories and
assigning the corresponding value to each key
to find the name, the ID from the attributes (9th) column must be split and assigned

the script finally writes the corresponding values appropriate for a BED format
'''
import argparse

def main():
    parser = argparse.ArgumentParser( description="Convert gff3 file to a bed file")
    
    parser.add_argument("-i", "--input_file", type=str, required=True, help="input file to be translated")
    parser.add_argument("-o", "--output_file", type=str, required=True, help="output file to be created")
    
    args = parser.parse_args()
    
    gff3File = open(args.input_file, "r", encoding="utf-8")
    bedFile = open(args.output_file, "w", encoding="utf-8")
    
    for line in gff3File:
        readLine = line.rstrip()
        if readLine.startswith("#"):
            continue
        else:
            elements = readLine.split("\t")
            columns = dict()
            columns["seqid"] = elements[0]
            columns["source"] = elements[1]
            columns["type"] = elements[2]
            columns["start"] = int(elements[3])
            columns["end"] = int(elements[4])
            columns["score"] = elements[5]
            columns["strand"] = elements[6]
            columns["phase"] = elements[7]
            columns["elements"] = elements[8]
            
            attributes = dict()
            attributes_list = columns["elements"].split(";")
            attributes_list2 = [x.split("=",1) for x in attributes_list if "=" in x]
            attributes = dict(attributes_list2)
            
            try:
                columns["name"] = attributes["ID"]
            except KeyError:
                columns["name"] = "."
                          
            columns["start"] -= 1

            if columns["type"] == "gene":
                columns["score"] = 0
                bedFile.write(columns["seqid"]+ "\t"+
                      str(columns["start"])+ "\t"+
                      str(columns["end"])+ "\t"+
                      columns["name"]+ "\t"+
                      str(columns["score"])+ "\t"+
                      columns["strand"]+
                              "\n")
    return 0
    gff3File.close()   
    bedFile.close()


if __name__ == "__main__":
    main()
