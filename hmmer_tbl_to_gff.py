import pandas #Para los dataframes
import argparse #Para el parseador de argumentos
import re

#Crear al parseador
parser = argparse.ArgumentParser() 

#Añadir argumentos
parser.add_argument("--file", type = str, required = True, help = "The input file produced by blast with output format 6 (tabular).")
parser.add_argument("--bitscore", type = float, required = True, help = "Bitscore obtained by blast. Used to add a level of filter to the result. If 0, all results will be maintained.")
parser.add_argument("--source", type = str, required = True, help = "The mode in which blast was performed, e.g. blastn, blastp, tblastx, ...")

#Parsear argumentos hacia una variable de almacenamiento
args = parser.parse_args()


#Renombrar variables o modificarlas si es necesario.
file = str(args.file)
outname = file.removesuffix('.tbl') + '.gff' #Cambiar nombre de salida
method = str(args.source)
bitscore = float(args.bitscore)


df = pandas.DataFrame(columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]) # Iniciar un df vacio pero con nombres de columnas

with open(file) as filehandler: #Abrir archivo
    for line in filehandler: #Para cada linea
        if line.startswith('#'):
            continue
        else:
            line = line.rstrip().replace("|", "; ") #Corta el trailing y remplaza | con ;
            line_elements = re.sub(r"\s+", " ", line).replace(' ', '\t').split('\t') #Los tbl tienen espacios en blanco. Se eliminan, se deja uno para reemplazar por tabsy se separa por tabs. Se crea una lista
            if float(line_elements[13]) > bitscore: #Si bitscore es mayor que el declarado en las opciones.
                new_line = line_elements[2] + "\t" + method + "\t" + "RSS" + "\t" + line_elements[6] + "\t" + line_elements[7] + "\t.\t.\t.\t" + "Query=" + line_elements[0] #Reordena la linea como gff. Se crea un string
                new_line = new_line.split("\t") #Partir el string por tabulacion
                if int(new_line[4]) > int(new_line[3]): #Si end es mayor que start
                    new_line[6] = "+" #Cambiar sexta posicion a +
                else: #Sino
                    new_line[4], new_line[3] = new_line[3], new_line[4] #Cambiar de posición los elementos de lista
                    new_line[6] = "-" #Cambiar sexto elemento a -
                df_list = pandas.DataFrame.from_dict({"seqid": [new_line[0]], "source": [new_line[1]], "type": [new_line[2]], "start": [new_line[3]], "end": [new_line[4]], \
                    "score": [new_line[5]], "strand": [new_line[6]], "phase": [new_line[7]], "attributes": [new_line[8]] } ) #Crea un df de un dicc, a partir de new line y asignando nombres de columnas respectivos.
                df = pandas.concat([df, df_list]) #Unir df que tienen los mismo nombres de columnas
df.to_csv(outname , sep = "\t", header = False,  index = False)
