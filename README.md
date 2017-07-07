### Как запускать

1) Подготовить окружение
```bash
$ docker build -t fenics/ice_island .

$ curl -s https://get.fenicsproject.org | bash
```

2) Далее запускать
```bash
$ chmod +x ./run.sh
$ ./run.sh --image=fenics/ice_island --config=path/to/config.json
```
