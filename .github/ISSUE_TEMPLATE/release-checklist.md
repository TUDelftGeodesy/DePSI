---
name: Release Checklist
about: Create a checklist for making a new release.
title: Release checklist vX.Y.Z
labels: ''
assignees: ''
type: Task

---

I will make a new release by preforming the following procedures:

- [ ]Step1: Update version in pyproject.toml

- [ ]Step2: Update `docs/CHANGELOG.md` with new change information

- [ ]Step3: Create release/tag at https://github.com/TUDelftGeodesy/sarxarray/releases

Hints:
- Step 1 and 2 can be done semi-automatically by calling agent skill in [release-changelog](.github/skills/release-changelog). However the changes will always be reviewed and commited by human.
