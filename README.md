# Nevroevolucija za strojno učenje iz tabelaričnih podatkov

Build programa lahko izvedete na naslednji način:
- ustvarite "release" CMake profil z ukazom `cmake -DCMAKE_BUILD_TYPE=Release -S/absolutna_pot/Neuroevolution -B/absolutna_pot/Neuroevolution/cmake-build-release`
- poženete build z ukazom `cmake --build /absolutna_pot/Neuroevolution/cmake-build-release --target Neuroevolution -- -j 4`
- če želite program razhroščevati, lahko ustvarite tudi "debug" profil 

Zbuildan program bo na lokaciji /absolutna_pot/Neuroevolution/cmake-build-release/Neuroevolution.

V datoteki CmakeLists.txt lahko spreminjate compiler zastavice.
Debug profil ima privzeto nastavljeno `-Og -ggdb3`, release pa `-O3`.
Oba profila imata nastavljeni zastavici `-march='native' -mtune='native'`,
zaradi katerih compiler čim bolje izkoristi procesorsko arhitekturo računalnika, na katerem poženemo build.
Zbuildan program pa ne bo deloval na procesorjih z različno arhitekturo.

Program ne potrebuje nobenih dodatno nameščenih knjižnic, dovolj je samo paket *libstdc++*.

Arhitekturo svojega računalnika in podprte optimizacijske možnosti lahko preverite z ukazom `gcc -march=native -mtune=native -Q --help=target`.

Program lahko poženete na naslednji način `./Neuroevolution car.data 350 40 100 4 true 0.1 175 true -0.00001 300 MCC`.
Več informacij o podprtih argumentih se nahaja v *main.cpp*.
Če želite ustvarjeno datoteko `graph.json` prikazati, jo lahko kopirate v mapo *visualization*, *visualization.html*
pa odprete v brskalniku. Datoteka `graphFull.json` vsebuje ne-minimiziran graf, prikažete pa jo lahko tako, da v 
*visualization.html* popravite vrstico 13.

Uporabljene podatkovne množice Iris, Wine, Car Evaluation in Statlog se nahajajo v mapi *datasets*, 
rezultati testiranja pa so zbrani v mapah *ods* in *orange*.
Izvorna koda poročila je v mapi *porocilo/proper*, za izdelavo PDF datoteke pa potrebujete nameščen program Inkscape,
s katerim se izvede rasterizacija vektorskih slik.