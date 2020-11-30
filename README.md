# Pigcasso

Inspired by the world-renowned Painting Pig &ndash; Pigcasso &ndash; is a Python utility plotting data from `.vtk` files (tipically the output of CFD simulations) to many other formats &ndash; such as `.eps` or `.png` &ndash; everything relies on `matplotlib`, thus for the moment just the possibility to plot `2D` and `1D` figures is available.

<img src="doc/gallery/pigcasso.jpg" alt="drawing" width="1000"/>

*Rescued from the slaughterhouse in May 2016, Pigcasso is one fine, fortunate swine* <!-- taken from here: https://www.vegangame.it/storie-di-animali/pigcasso-la-maialina-dallestro-creativo -->

At almost 2000 pounds, it's the heavy-weight Abstarct Expressionist of the World! Have a look at its [masterpieces](https://pigcasso.myshopify.com).

## How to make pigcasso paint for you?

The repo workflow is organised as follows:

| branch        | |
| :-------------: |-------------| 
| `master` | where it goes the code validated after the production phase | 
| `dev` | where it goes the code deployed after testing ready for production | 
| `backport` | where it goes the code backported from other container repos | 

To use `pigcasso` we recommend to install it "cloning" its content within a `container` repo &ndash; the one providing for the data to be plotted &ndash; and the best way to do this is by merging it as a `subtree`. Have a look at [this guide](https://medium.com/@porteneuve/mastering-git-subtrees-943d29a798ec) (and the following `git 2.9` updates within the comments, there below) if not refer to the GitHub [official guide](https://help.github.com/en/github/using-git/about-git-subtree-merges).

## How to get updates from pigcasso?

Assuming you are getting an update from the subtree’s remote, according to [the guide](https://medium.com/@porteneuve/mastering-git-subtrees-943d29a798ec), after checking out your feature branch, on the `container` repo, most of the times `git merge -s subtree --allow-unrelated-histories --squash pigcasso/backport` will be enough, if not opt for the default strategy (recursive) with an explicit prefix through its subtree option:
```
git merge -X subtree=3rd-party/pigcasso/ --allow-unrelated-histories --squash pigcasso/dev
``` 
and eventually solve any conflict that may arise. It can happen that the heuristics used by the subtree merge strategy to figure out the subdirectory prefix get confused.

## How to backport to pigcasso?

To backport it's quite an easy task, according to [the guide](https://medium.com/@porteneuve/mastering-git-subtrees-943d29a798ec), after pushing "the" commit (to backport here) on the `container` repo, you will have to follow these steps:
```
git stash ; wrk_branch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
git remote add -f pigcasso https://github.com/andreagalle/pigcasso.git
git checkout -b oink pigcasso/backport
git pull pigcasso backport
git cherry-pick -x --strategy=subtree $wrk_branch && git show HEAD --stat
git push pigcasso oink:backport
git checkout $wrk_branch && git stash apply
```
of course the `git show HEAD --stat` command is not neccessary, but advised to check everything before pushing upstream here. To check in more detail `git log --graph --decorate --oneline` and `git show HEAD` are always the best options.

