# Compilateur utilisé
CC=mpicxx

# Options en mode debug - La variable est DEBUG est définie comme vraie
DEBUG_FLAG = --std=c++11

# On choisit comment on compile
CXX_FLAGS = $(DEBUG_FLAG)

# Le nom de l'exécutable
PROG = exec

# Les fichiers source à compiler
SRC = main.cc TimeScheme.cpp FiniteVolume.cpp save.cpp function.cpp DataFile.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CC) $(CXX_FLAGS) -o $(PROG) $(SRC)
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG)
