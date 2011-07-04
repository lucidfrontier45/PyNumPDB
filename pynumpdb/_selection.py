from _sequence import codes_three

"""
A package for selecting atoms in PDB format

PDB format

 COLUMNS     DATA TYPE        FIELD         DEFINITION
 --------------------------------------------------------------
   1 - 6      Record name      "HETATM"
   7 - 11     Integer          serial        Atom serial number.
  13 - 16     Atom             name          Atom name.
  17          Character        altLoc        Alternate location indicator.
  18 - 20     Residue name     resName       Residue name.
  22          Character        chainID       Chain identifier.
  23 - 26     Integer          resSeq        Residue sequence number.
  27          AChar            iCode         Code for insertion of residues.
  31 - 38     Real(8.3)        x             Orthogonal coordinates for X.
  39 - 46     Real(8.3)        y             Orthogonal coordinates for Y.
  47 - 54     Real(8.3)        z             Orthogonal coordinates for Z.
  55 - 60     Real(6.2)        occupancy     Occupancy.
  61 - 66     Real(6.2)        tempFactor    Temperature factor.
  77 - 78     LString(2)       element       Element symbol; right-justified.
  79 - 80     LString(2)       charge        Charge on the atom.
"""

heavy_atoms = ("N","CA","C","O","CB","OG","SG","CG1","OG1","CG2","CD","SD",\
    "OD1","CD1","ND1","CE","OD2","ND2","CD2","NE")

def stripHydrogen(pdb_data):
  new_data = []
  for line in pdb_data:
    rec_name = line[0:6].strip()
    if rec_name in ("ATOM","HETATM"):
      atom_name = line[12:16].strip()
      if (("H" not in atom_name) or (atom_name in ("NH1","NH2","OH","CH2"))):
        new_data.append(line)
    else :
      new_data.append(line)
  return new_data

def getLoc(pdb_data,loc="A"):
  new_data = []
  for line in pdb_data:
    rec_name = line[0:6].strip()
    if rec_name in ("ATOM","HETATM"):
      altLoc = line[16]
      if altLoc in (" ",loc) :
        new_line = line[:16] + " " + line[17:]
        new_data.append(new_line)
    else :
      new_data.append(line)
  return new_data

def getModel(pdb_data,model_num=1,noh=True,loc=None):
  if noh:
    data = stripHydrogen(pdb_data)
  else:
    data = pdb_data
  
  if loc:
    data = getLoc(data,loc)

  new_data = []
  model_count = 0
  #first_chain = None

  for line in data:

    # choose first model
    rec_name = line[0:6].strip()
    if rec_name == "MODEL" :
      model_count += 1
    if model_count < model_num :
      continue
    elif model_count > model_num :
      break

    if rec_name not in ("ATOM","HETATM"):
      continue
      
    #if first_chain == None:
    #  first_chain = line[21]
    #if line[21] != first_chain:
    #  break

    new_data.append(line)

  if model_count == 0:
    new_data = data

  return new_data

def _getChain(pdb_data,ids=("A",)):

  no_id = True
  new_data = []
  for line in pdb_data:

    rec_name = line[0:6].strip()
    if rec_name not in ("ATOM","HETATM"):
      continue

    chain_id = line[21]
    if chain_id != " ":
      no_id = False
    if chain_id in ids:
      new_data.append(line)

  if no_id:
    return pdb_data,(" ",)
  else:
    return new_data,ids

def _getFirstChain(pdb_data):
  chain_id = ""
  for line in pdb_data:
    rec_name = line[0:6].strip()
    if rec_name in ("ATOM","HETATM"):
      res_name = line[17:20]
      if res_name in codes_three:
        chain_id = line[21]
        break
  return _getChain(pdb_data,chain_id)

def getChain(pdb_data,ids=None):
  if ids == None:
    return _getFirstChain(pdb_data)
  else:
    return _getChain(pdb_data,ids)

def reNumRes(pdb_data,start=1):
  iatom = 0
  lastAtom = "dummy"
  lastResId = -999
  
  ires = int(start) - 1
  new_data = []
  
  for line in getLoc(pdb_data):
   
    recName = line[0:6].strip()
    if not recName in ("ATOM","HETATM"):
      continue 

    resId = int(line[22:26].strip())
    if not resId == lastResId:
      lastResId = resId
      ires += 1
  
    iatom += 1
    new_data.append("%s%5d%s%4d%s" % (line[0:6],iatom,line[11:22],ires,line[26:]))

  return new_data

def selectRes(pdb_data,start,end):
  new_data = []
  for line in pdb_data:

    recName = line[0:6].strip()
    if not recName in ("ATOM","HETATM"):
      continue

    res_id = int(line[22:26].strip())
    if res_id < start:
      continue
    elif res_id > end:
      break
    else:
      new_data.append(line)

  return new_data

def getPeptide(pdb_data,start,end,reReNum=False,noh=True):
  model,ids = getChain(getModel(pdb_data,noh=noh))
  selected = selectRes(reNumRes(model),start,end)
  if reReNum:
    return reNumRes(selected)
  else :
    return selected

