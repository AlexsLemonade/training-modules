# An introductory guide to debugging 

## Tips for approaching error

### 1) Identify which phrase is the source of the problem 

### 2) Look at the documentation for a function to make sure you are using it correctly

#### Use the help bar

#### If you are using a Bioconductor package, look at their documents

#### a) What structure object is the function built to take? 

#### b) What character type is the function built to take?  

### 3) Google it

#### StackOverflow
#### Bioconductor 
#### GitHub Issues

### 4) Google it again

### 5) Look at the source code for that function

## A guide to the most common errors

```
Error in file(file, "rt") : cannot open the connection
In addition: Warning message:
In file(file, "rt") : cannot open file '<FILENAME>': No such file or directory
```
```
Error in ... could not find function <FUNCTION_NAME>
```

```
Error in ... object '<OBJECT_NAME>' not found
```

```
Error in library(<PACKAGE_NAME>) : there is no package called ‘<PACKAGE_NAME>’
```

```
Error in ... no applicable method for <OBJECT_NAME> applied to an object of class <CLASS_OF_OBJECT>
```

```
Error in ... subscript out of bounds
```

```
ERROR: compilation failed for package <PACKAGE_NAME>
```
