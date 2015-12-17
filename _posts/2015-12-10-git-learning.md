---
layout: post
title:  "git learning"
date:   2015-12-10
categories: git-workflow
---

Here are some notes (mostly for myself) about common workflows when using [git][git-link], from what I have learned on [Codecademy][codecademy-git] lesson on git. I am now frequently using [GitHub][github] to work on my projects. I am just getting into branching, and hopefully referring to this will help with simplifying the process of testing and developing code/scripts, and help when collaborating with others.

```
	## standard workflow
	git add filename
	git commit -m "message"
	git show HEAD # this shows the most recent commit
	git checkout HEAD filename # This restores the file to last commit
	git reset HEAD filename # This unstages a commit
	git reset ####### use the first 7 characters of the commit SHA in the Git log

	## branches
	git branch # check the branch you are on
	git branch new_branch # add a branch
	git checkout new_branch # switch to new branch
	git branch # check branch again
	
	## edits with branches
	git checkout master # switch to master branch
	git merge new_branch # this will `fast-forward' the master to the new branch
	# if conflicts arise: edit the file and remove git markings
	git branch -d branch_to_be_deleted

	git clone remote_location local_name
	git remote -v # list the remotes, one for fetch, one for push
	git fetch remote_location # this DOES NOT MERGE -- it makes a remote branch
	git merge origin/master # merge/fast-forward the origin/master remote branch to local
	
	#1. git clone remote_location local_name
	#2. git branch local_branch_name
	#3. git add file; git commit -m "message" # develop feature, commit
	#4. git fetch remote_location; git merge origin/master
	#5. git push origin local_branch_name
	
```

In words, the workflow with branches is the following:

1. fetch and merge from remote
2. create a branch
3. develop feature, commit
4. fetch and merge from remote
5. push branch up to remote for review

Steps 1 and 4 protect against merge conflicts. If you keep up to date and no one else is editing your code, you may skip these steps.

[git-link]: https://git-scm.com
[codecademy-git]: https://www.codecademy.com/learn/learn-git
[github]: https://training.github.com/classes/introduction/