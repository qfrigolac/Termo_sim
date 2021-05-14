# Simulació - Termodinàmica i Mecànica Estadística
<img align="left" src="https://img.shields.io/badge/Termodinàmica-Simulació-yellow"> <img align="left" src="https://img.shields.io/badge/Llenguatge-Python-blue"> <br>


Codi de Bruce Sherwood. Modificat per Joaquim Frigola i Pau Sànchez.

## Introducció
En aquest codi se simula unes partícules d'un gas per tal d'extreure'n diverses conclusions.

## Variables
Al principi del codi trobem trobem algunes variables que podem modificar per fer les diferents simulacions.

| Variable | Funció                                                                                                               |
| :----- | :------------------------------------------------------------------------------------------------------------------- |
| Natoms   | Nombre d'àtoms (boles) a la simulació                                                                                |
| L        | Mida de la caixa (lateral)                                                                                           |
| Ratom    | Radi dels àtoms                                                                                                      |
| k        | Constant de Boltzmann                                                                                                |
| T0       | Temperatura en el qual es realitzarà la simulació                                                                    |
| dt       | Quantitat de temps que s'avança en cada pas                                                                          |
| nderiv   | Nombre de iteracions que farà el programa abans de registrar les magnituds termodinàmiques al calcular les derivades |
| deltaL   | Augment de la caixa en fer les derivades                                                                             |
| dq4      | Augment de la caixa en l'expansió lliure  en l'eix x                                                                 |
## Funcionament
Al prémer el botó "Fer derivades programades" el programa calcularà 28 simulacions per computar les derivades i el coeficient de dilatació. Podem parar-ho apretant el mateix botó. Si premem "Restart" es recomençarà la simulació amb els valors per defecte.

Al prémer "Expansió lliure (pas)" es reinicialitzarà la simulació amb la caixa modificada en l'eix x, amb un increment dq4. El botó "Instantània pV" guarda un punt de pressió i invers del volum i ho representa a un gràfic.

Finalment, "R. Pressió Mitjana" esborra la mitjana global de la pressió (línia carbassa) i la torna a començar des de zero.

## Informació disponible
Quan fem la simulació a la pantalla se'ns mostra informació diversa:

| Informació          | Significat                                                                                                                                                                                                           |
| :------------------ | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Pressió mitjana     | Mitjana de les pressió calculades segons les partícules que xoquen a les parets. S'actualitza cada múltiple de 10 d'n.                                                                                               |
| Temperatura mitjana | Temperatura calculada amb la mitjana de les velocitats de cada partícula                                                                                                                                             |
| pV                  | Pressió mitjana per volum                                                                                                                                                                                            |
| NkT                 | Nombre total d'àtoms per la constant de Boltzmann per la temperatura mitjana                                                                                                                                         |
| n                   | Nombre de passos des de l'inici                                                                                                                                                                                      |
| Gràfic 1            | Histograma del nombre d'àtoms en funció de velocitat (en vermell) i distribució teòrica en l'equilibri (blau)                                                                                                        |
| Gràfic 2            | Pressió instantània en funció del temps (en blau), mitjana local d'aquesta pressió per cada n=10 (en groc) i mitjana global de la pressió des de l'inici de la simulació o des de "R. Pressió Mitjana" (en carbassa) |
| Gràfic 3            | Temperaturatura calculada amb vx en funció del temps (en groc) i temperatura calculada usant la mitjana de les velocitats (en blau)                                                                                  |
| Gràfic 4            | Histograma del coeficient de dilatació calculat                                                                                                                                                                      |
| Gràfic 5     |   Pressió en funció de l'invers del volum quan es prem "Instantània pV"                                                                                                                |

Noteu que els dos últims gràfics poden sortir girats.
