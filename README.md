# ABModelFunction

Fitness function used to evaluate AB model that represents proten.
Solution vector is defined with two angles.

## Example of the input vector
 
If the input sequence is ABABABA:

* Input sequence length is 7
* Input vector `x` dimensionalty is equal to 10
* First 5 varibles of vector `x` represent first angle
* Last 5 varibles of vector `x` represent the second angle

## Example of usage

### Using predefined proterin in the package database

```python
import numpy as np
from abmodel.abmodel import Model
from abmodel.abmodel import pfo_instance

# Get a known sequence from database
seq, opt = pfo_instance('2EWH')
print(seq, ' = ', opt)

# Generate the model function
m = Model(seq)

# Prepare date for processing
lb, ub = m.bounds()
lb, ub = np.asarray(lb), np.asarray(ub)
dim = m.dim()

# Generate a random solution x and evaluate it
x = lb + np.random.rand(dim) * (ub - lb)
print(m.eval(x))
```
### Using user define AB model

```python
import numpy as np
from abmodel.abmodel import Model

# User defined sequence
seq = 'ABBBABBABBBABBBAAAB'

# Generate the model function
m = Model(seq)

# Prepare date for processing
lb, ub = m.bounds()
lb, ub = np.asarray(lb), np.asarray(ub)
dim = m.dim()

# Generate a random solution x and evaluate it
x = lb + np.random.rand(dim) * (ub - lb)
print(m.eval(x))
```