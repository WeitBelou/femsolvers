### Как запускать

1) Подготовить окружение
```bash
$ docker build -t fenics/ice_island .

$ curl -s https://get.fenicsproject.org | bash
$ fenicsproject create ice_island fenics/ice_island
```

2) Далее запускать
```bash
$ fenicsproject start ice_island
$ python3 src/runner.py path/to/config.yml
```
