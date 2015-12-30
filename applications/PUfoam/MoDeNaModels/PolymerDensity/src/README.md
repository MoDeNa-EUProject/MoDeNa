# Readme file for PCSAFTi\_Density

The program is automatically compiled and linked by executing the following two commands:

```
cmake .
make
```

Afterwards it can be used in some way that Jonas will now explain

# Minimal test to check that the program works

Create a file with inputs to the program:

```
printf "270.0\n2\nair\npu\n0.\n0." > in.txt
```

Create a file in which the program can save the output

```
touch out.txt
```

```
./PCSAFT_Density
```
