# multi agents simulation
# Exercice 2
## Exercice 2_1
```julia
helium = Molecule([0.0,0.0,0.0],[0.0,0.0,0.0],6.6464731e-27,3.1e-11,"He",[],[])
Neon = Molecule([0.0,0.0,0.0],[0.0,0.0,0.0],3.3509177e-26,4.8e-11,"Ne",[],[])
diazote = Molecule([0.0,0.0,0.0],[0.0,0.0,0.0],1.85e-10,4.65e-26,"N2",[],[])
dioxygene = Molecule([0.0,0.0,0.0],[0.0,0.0,0.0],1.5e-10,5.3135e-26,"O2",[],[])
```
## Exercice 2_2
Un gaz 

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
