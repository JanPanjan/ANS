![[Pasted image 20250424161130.png]]

>[!bug] DODAJ V GITIGNORE VSE VRSTE DATOTEK V ANS IN UVB DIRECTORY

![[Pasted image 20250424161230.png]]

# 1. dodaj v .gitignore

Če ga še nisi ustvaril, ga ustvari in dodaj not vse te ogromne datoteke in nato ga dodaj v commit.

```bash
git add .gitignore
```

Če ga že imaš, dodaj imena datotek.

# 2. odstrani iz commit-a

```bash
git rm --cached -r <dir>
```
 
- `git rm --cached` odstraani datoteke iz commita
- `-r` pomeni recursive

# 2. update (ammend-aj) commit

```bash
git commit --amend -m "something (excluding large files)"`
```

- `git commit --amend` omogoča da popraviš zadnji commit
- z `-m` posodobiš commit message

# 3. push

Probaj potisnit spremembe.

```bash
git push
```