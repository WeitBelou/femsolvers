### Как запускать

1) Подготовить окружение
```bash
$ curl -s https://get.fenicsproject.org | bash
$ fenicsproject create ice_island

$ fenicsproject start ice_island
$ pip3 install -r requirements.txt
$ exit
```

2) Далее запускать
```bash
$ fenicsproject start ice_island
$ python3 src/runner.py path/to/config.yml
```
