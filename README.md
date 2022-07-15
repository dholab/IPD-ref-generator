# IPD-ref-generator

```
chmod +x bin/*.py
```

```
docker build --tag ipd-ref-generator:v1_0_5 config/
docker run -it -v $(pwd)/:/scratch -w /scratch ipd-ref-generator:v1_0_5
nextflow run ipd-ref-generator.nf
```

Note that to build the above docker container, you may need to increase the amount of memory alloted to Docker in the Docker Engine preferences.
