#On purge la liste des suffixes utilis� pour les r�les implicites
.SUFFIXES:

#On ajoute simplements les extensions dont l'on a besoin
.SUFFIXES:.cpp .o

#Nom de l'executable
EXEC=tp6

#Liste des fichiers sources separes par des espaces
SOURCES=main.cpp

#Liste des fichiers objets
OBJETS=$(SOURCES:%.cpp=%.o)

#Compilateur et options de compilation
CCPP=g++
CFLAGS= -W -Wall -Wextra -Werror -Wno-unused-parameter -pedantic -std=c++11
LFLAGS= -L . -lm -lgmp

#R�le explicite de construction de l'ex�utable
$(EXEC):$(OBJETS) Makefile
	$(CCPP) -o  $(EXEC) $(OBJETS) $(LFLAGS)
.cpp.o:
	$(CCPP) $(CFLAGS) -c $< -o $@

clean:
	rm $(OBJETS)
clear:
	rm $(EXEC)
depend:
	sed -e "/^#DEPENDANCIES/,$$ d" Makefile >dependances
	echo "#DEPENDANCIES" >> dependances
	$(CCPP) -MM $(SOURCES) >> dependances
	cat dependances >Makefile
	rm dependances

#DEPENDANCIES
main.o: main.cpp
