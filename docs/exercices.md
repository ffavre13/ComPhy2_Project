# multi agents simulation
# Exercice 2
## Exercice 2_1
### Masse
Pour calculer la masse, on utilise la masse molaire pour chaque gaz et on la divise par $N_A$ (nombre d’Avogadro) pour obtenir la masse d’une molécule en [kg]. Le 1000 c'est pour passer de [g] en [kg]. Si on a N2, faire la masse molaire présente dans le tableau périodique * 2.

$$
m = \frac{M}{N_A * 1000}
$$

- M = masse molaire [g/mol]
- $N_A$ = $6.022*10^{23}$ [1/mol]

- He
    - Masse molaire = 4.0026
    - Masse [kg] = 6.64663e-27
- Ne
    - Masse molaire = 20.180
    - Masse [kg] = 3.35105e-26
- N2
    - Masse molaire = 14.007 * 2 = 28.014
    - Masse [kg] = 4.65194e-26
- O2
    - Masse molaire = 15.999 * 2 = 31.998
    - Masse [kg] = 5.31352e-26

### Rayon

On prends le rayon de van der Waals : https://en.wikipedia.org/wiki/Van_der_Waals_radius

- He
    - https://fr.wikipedia.org/wiki/H%C3%A9lium
    - Rayon atomique [m] : 1.4e-10 
- Ne
    - https://fr.wikipedia.org/wiki/N%C3%A9on
    - Rayon atomique [m] : 1.54e-10 
- N2
    - https://en.wikipedia.org/wiki/Van_der_Waals_radius tableau a la fin de la page
    - ~ le rayon de N
    - Rayon atomique [m] : 1.55e-10
- O2
    - https://en.wikipedia.org/wiki/Van_der_Waals_radius tableau a la fin de la page
    - ~ le rayon de O
    - Rayon atomique [m] : 1.52e-10 
### Formule chimique

- He
- Ne
- N2
- O2

```julia
helium = Molecule([0.0,0.0,0.0],[0.0,0.0,0.0],6.64663e-27,1.4e-10 ,"He",[],[])
Neon = Molecule([0.0,0.0,0.0],[0.0,0.0,0.0],3.35105e-26,1.54e-10 ,"Ne",[],[])
diazote = Molecule([0.0,0.0,0.0],[0.0,0.0,0.0],4.65194e-26,1.55e-10,"N2",[],[])
dioxygene = Molecule([0.0,0.0,0.0],[0.0,0.0,0.0],5.31352e-26,1.52e-10 ,"O2",[],[])
```
## Exercice 2_2
La différentes est que on a He et Ne qui sont des atomes et N2 et O2 qui sont des molécules composée de 2 atomes. Comme la liason covalente dans les molécules N2 et O2 est très petite donc les atomes sont très proche, on peut considérer ces 2 molécules comme une sphère pour notre modèle. Donc la théorie s'applique bien avec cette lègre aproximation qui est de considérer que ça forme une sphère.

# Exercice 3

## Exercice 3_1

Il n'y a pas de force entre les molécules, donc il n'y a pas d'accélération. Les seuls intéractions entre les molécules sera lors des collisions

$
m \cdot \overrightarrow{a} = 0
$

$
\overrightarrow{a} = 0
$

$
\overrightarrow{v} = \overrightarrow{v_0}
$

$
\overrightarrow{r} = \overrightarrow{r_0} + \overrightarrow{v} (t)
$

## Exercice 3_2

$
\overrightarrow{r}(t + \Delta t) = \overrightarrow{r}(t) + \Delta t \cdot \overrightarrow{v} (t) 
$

# Exercice 4

## Exercice 4_1

On a le modèle des gaz parfais donc les molécules molécules sont assimilées à des sphères. Elles n’interagissent pas à distance.
Lorsqu’elles se touchent, elles subissent un choc élastique.

## Exercice 4_2

### Scenarios
#### 1

Les 2 molécules se percutent avec la même masse, la vitesse opposée, et la direction opposée.
Normalement, elle devrait se repousser

```julia
positions = [[-0.1,0.0,0.0],[0.1,0.0,0.0]]
velocities = [[5.0,0.0,0.0],[-5.0,0.0,0.0]]
masses = [1.0,1.0]
radius = [0.01,0.01]
chimical_formulas = ["TEST","TEST"]
```

#### 2

Il y a une molécule qui est immobile et une qui lui fonce dessus. Elles ont la même masses.
Normalement, la molécule immobile devrait prendre toute la vitesse de cela qui la percute et cela qui était en mouvement se stoppe. 

```julia
positions = [[0.0,0.0,0.0],[0.1,0.0,0.0]]
velocities = [[0.0,0.0,0.0],[-5.0,0.0,0.0]]
masses = [1.0,2.0]
radius = [0.01,0.01]
chimical_formulas = ["TEST","TEST"]
```

#### 3

Les 2 molécules se percutent avec la même masse, la vitesse opposée, et la direction opposée avec un légé décalage en hauteur (rayon//2).
Normalement, elle devrait se repousser et partir avec une angle de 45°

```julia
positions = [[-0.1,0.0,0.005],[0.1,0.0,0.0]]
velocities = [[5.0,0.0,0.0],[-5.0,0.0,0.0]]
masses = [1.0,1.0]
radius = [0.01,0.01]
chimical_formulas = ["TEST","TEST"]
```

#### 4

Les 2 molécules se percutent avec la une masse différente, la vitesse dans le même sens mais pas la même, celle devant a une vitesse supérieur, et la direction opposée avec un légé décalage en hauteur (rayon//2).
Normalement, elle devrait se repousser et partir avec une angle de 45°

```julia
positions = [[0.0,0.0,0.0],[0.1,0.0,0.0]]
velocities = [[-2.5,0.0,0.0],[-10.0,0.0,0.0]]
masses = [0.000001,1.0]
radius = [0.01,0.01]
chimical_formulas = ["TEST","TEST"]
```
