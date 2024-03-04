# AtmoTerre

Ce code permet de calculer la spéciation d'une atmosphère en équilibre avec un océan magmatique.

## Utilisation

Pour utiliser ce code, exécutez le script `atmoter.py` avec l'une des commandes suivantes :

- Pour imprimer les données des tables JANAF à une température donnée

```bash
python3 atmoter.py --t 298.15
```

- Pour calculer la fugacité d'oxygène du magma.

```bash
python3 atmoter.py fO2
```

- Pour calculer la spéciation de l'atmosphère.

```bash
python3 atmoter.py speciation
```

- Pour calculer une courbe de spéciation.

```bash
python3 atmoter.py curve
```

### Options

Vous pouvez également spécifier les options suivantes :

- `--t`
    Fixe la température en Kelvin.
- `--p`
    Fixe la pression en bar.
- `--ho`
    Fixe le rapport atomique H/O.

## Auteur

- Régis THIERY

## Contact

Pour toute question ou suggestion, veuillez contacter Régis THIERY à l'adresse e-mail suivante : <regthiery@gmail.com>
