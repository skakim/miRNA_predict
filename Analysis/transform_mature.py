with open("mature.fa","r") as mature_input:
    with open("filtered_mature.fa","w") as mature_output:
        hsa = False
        for line in mature_input:
            if line[0] == ">":
                if "hsa" in line:
                    hsa = True
                    mature_output.write(line[1:].replace("\n"," "))
            elif hsa:
                mature_output.write(line)
                hsa = False
    mature_output.close()
mature_input.close()
