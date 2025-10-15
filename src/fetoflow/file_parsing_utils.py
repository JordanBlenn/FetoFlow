def read_nodes(filename):
    if not isinstance(filename,str):
        raise(TypeError("File name should be a string!"))
    if not filename[-7:] == ".ipnode":
        raise(TypeError("This function expects a .ipnode file."))
    with open(filename, "r") as f:
        lines = f.readlines()
        nodes = {}
        numVars = int(lines[4].split()[-1])
        baseStep = numVars + 2
        i = 4 + 2 * numVars + 2
        while i < len(lines):
            node_id = int(lines[i].split()[-1])
            nversions = int(lines[i + 1].split()[-1])
            i_step = baseStep + nversions * numVars + (nversions * numVars if nversions > 1 else 0)
            c_step = numVars * nversions - 1
            coords = []
            for c in range(1 + nversions, i_step, c_step):
                coord = lines[i + c].split()[-1]
                coords.append(float(coord))
            nodes[node_id - 1] = coords  # 0-based indexing for the networkX geometry.
            i += i_step
    return nodes


def read_elements(filename):
    if not isinstance(filename,str):
        raise(TypeError("File name should be a string!"))
    if not filename[-7:] == ".ipelem":
        raise(TypeError("This function expects a .ipelem file."))
    
    with open(filename, "r") as f:
        lines = f.readlines()
        elems = []
        i = 5
        while i < len(lines):
            intraElementStep = 5
            nodes = tuple(int(x) - 1 for x in lines[i + intraElementStep].split()[-2:])  # Translate to 0-based indexing.
            elems.append(nodes)
            while len(lines[i].split()) != 0:
                i += 1
                if i + 1 == len(lines):
                    return elems
            i += 1
        return elems


def define_fields_from_files(files: dict[str]):
    """
    Defines field(s) as specified in ipfield file(s).
    """
    if not isinstance(files, dict):
        raise (TypeError("files must be a dictionary in the format files[field_name] = filename"))
    fields = {}
    for field in files.keys():
        file_name = files[field]
        if not file_name[-7:] == ".ipfiel":
            ext_start = -(str.__reversed__(file_name).find(".") + 1)
            if ext_start is not None:
                raise(TypeError(f"This function expects a .ipfiel file, got {file_name[ext_start:]}"))
            else:
                raise(TypeError(f"This function expects a .ipfiel file. No file extension found."))
        with open(file_name, "r") as f:
            i = 7  # ignore metadata
            lines = f.readlines()
            max_digits = len(lines[3][lines[3].find(":") + 1 :].strip())
            currentField = {}
            while i < len(lines):
                j = i + 2
                id = (
                    int(lines[i][-max_digits-1:].strip()) - 1
                )  # assuming for now these correspond to element ids, -1 to 0 based
                val = float(lines[j][lines[j].find(":")+1:].strip())
                currentField[id] = val
                i += 4
        fields[field] = currentField
    return fields
