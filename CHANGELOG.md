# Changelog

## [1.6.0](https://github.com/xsitarcik/wrappers/compare/v1.5.8...v1.6.0) (2023-06-30)


### Features

* added basic samtools rules, faidx, index, sort and view ([a188213](https://github.com/xsitarcik/wrappers/commit/a1882135ea40ca48e9c55cdf605477a3a3543a27))
* added bowtie2 mapping wrappers ([cde91e4](https://github.com/xsitarcik/wrappers/commit/cde91e403af67ae05f60b3ded0b61c0a6413853f))
* added bwa mapping wrapper rules ([58e24d1](https://github.com/xsitarcik/wrappers/commit/58e24d18be0dd26e774829b11b3b4fb84c37f786))
* added cutadapt for both single and paired ([0ce8037](https://github.com/xsitarcik/wrappers/commit/0ce803703acfaad3b45486b5e066019120ea8bdd))
* added fastqc quality report ([878ca42](https://github.com/xsitarcik/wrappers/commit/878ca425aa9d2d32d09eb1df03815ff7c66c781e))
* added picard rules for BedToIntervalList, MarkDuplicates and CreateSequenceDictionary ([dd384ea](https://github.com/xsitarcik/wrappers/commit/dd384eac4649dd88f10f57af745867bacd3a5ea2))
* added qualimap bamqc basic rule ([87b053c](https://github.com/xsitarcik/wrappers/commit/87b053ca1fee3e956edfcf9d7f14edea14e6911a))
* added subsampling using seqtk ([a93ef17](https://github.com/xsitarcik/wrappers/commit/a93ef177fed8f1e1abcf427fa7f7cc9c8fe7928b))


### Bug Fixes

* added bwa_map with piped sort and optional filter ([#27](https://github.com/xsitarcik/wrappers/issues/27)) ([4f6071f](https://github.com/xsitarcik/wrappers/commit/4f6071f5c10dbed88710784dcfec87db9a980810))
* added temp dir for fastqc io java ([38a7d2f](https://github.com/xsitarcik/wrappers/commit/38a7d2fbf83e2e3ad3b18847463e168f13402875))
* added tmpdir handling in fastqc calling ([108768a](https://github.com/xsitarcik/wrappers/commit/108768af7e47960cb3c066e4dcce3a1a60131a34))
* bowtie2 wrapper to use correct wrapper repo instead of old one ([81bbc56](https://github.com/xsitarcik/wrappers/commit/81bbc56ed847756a1699e9be095c97b4fd82d993))
* bwa wrappers fixed to not use wildcards ([05ae11f](https://github.com/xsitarcik/wrappers/commit/05ae11f04f013aca0a14badc3d13e1ed07fced3b))
* bwa_mapping wrapper to use correct wrapper repo ([577142b](https://github.com/xsitarcik/wrappers/commit/577142b1fc345f0d996052a59408bcbf6dc22f78))
* cutadapt optional params removed and replaced by extra ([752cac6](https://github.com/xsitarcik/wrappers/commit/752cac6db1e143292313bbea255d79960d42ce5b))
* downgraded picard env ([5456d40](https://github.com/xsitarcik/wrappers/commit/5456d40ba5486755bfb83270c32085f125cca6cd))
* fastqc can be called with the memory parameter ([f79fb4a](https://github.com/xsitarcik/wrappers/commit/f79fb4aeeefe9c340d05c9e972c66ad049faf9e4))
* fastqc now correctly computes needed memory ([18d7479](https://github.com/xsitarcik/wrappers/commit/18d7479dc4a295e3fd5b642dc89744a59a47945f))
* fix seqtk memory mode argument ([d6ccb5b](https://github.com/xsitarcik/wrappers/commit/d6ccb5b9b087d6597acef1372dbb25a0d792d5e9))
* move to new repo ([98925b2](https://github.com/xsitarcik/wrappers/commit/98925b23aa5042c0c26176623a15f99f0cf61ee1))
* picard markduplicates using java mem arg ([f8b45cd](https://github.com/xsitarcik/wrappers/commit/f8b45cd0837f0c849e442a4d789d04e08f42bae1))

## [1.5.8](https://github.com/xsitarcik/wrappers/compare/v1.5.7...v1.5.8) (2023-06-27)


### Bug Fixes

* picard markduplicates using java mem arg ([f8b45cd](https://github.com/xsitarcik/wrappers/commit/f8b45cd0837f0c849e442a4d789d04e08f42bae1))

## [1.5.7](https://github.com/xsitarcik/wrappers/compare/v1.5.6...v1.5.7) (2023-06-14)


### Bug Fixes

* added bwa_map with piped sort and optional filter ([#27](https://github.com/xsitarcik/wrappers/issues/27)) ([4f6071f](https://github.com/xsitarcik/wrappers/commit/4f6071f5c10dbed88710784dcfec87db9a980810))

## [1.5.6](https://github.com/xsitarcik/wrappers/compare/v1.5.5...v1.5.6) (2023-06-13)


### Bug Fixes

* fastqc now correctly computes needed memory ([18d7479](https://github.com/xsitarcik/wrappers/commit/18d7479dc4a295e3fd5b642dc89744a59a47945f))

## [1.5.5](https://github.com/xsitarcik/wrappers/compare/v1.5.4...v1.5.5) (2023-06-13)


### Bug Fixes

* added temp dir for fastqc io java ([38a7d2f](https://github.com/xsitarcik/wrappers/commit/38a7d2fbf83e2e3ad3b18847463e168f13402875))

## [1.5.4](https://github.com/xsitarcik/wrappers/compare/v1.5.3...v1.5.4) (2023-06-09)


### Bug Fixes

* fastqc can be called with the memory parameter ([f79fb4a](https://github.com/xsitarcik/wrappers/commit/f79fb4aeeefe9c340d05c9e972c66ad049faf9e4))

## [1.5.3](https://github.com/xsitarcik/wrappers/compare/v1.5.2...v1.5.3) (2023-05-19)


### Bug Fixes

* added tmpdir handling in fastqc calling ([108768a](https://github.com/xsitarcik/wrappers/commit/108768af7e47960cb3c066e4dcce3a1a60131a34))

## [1.5.2](https://github.com/xsitarcik/wrappers/compare/v1.5.1...v1.5.2) (2023-05-11)


### Bug Fixes

* fix seqtk memory mode argument ([d6ccb5b](https://github.com/xsitarcik/wrappers/commit/d6ccb5b9b087d6597acef1372dbb25a0d792d5e9))

## [1.5.1](https://github.com/xsitarcik/wrappers/compare/v1.5.0...v1.5.1) (2023-05-11)


### Bug Fixes

* cutadapt optional params removed and replaced by extra ([752cac6](https://github.com/xsitarcik/wrappers/commit/752cac6db1e143292313bbea255d79960d42ce5b))

## [1.5.0](https://github.com/xsitarcik/wrappers/compare/v1.4.0...v1.5.0) (2023-05-10)


### Features

* added subsampling using seqtk ([a93ef17](https://github.com/xsitarcik/wrappers/commit/a93ef177fed8f1e1abcf427fa7f7cc9c8fe7928b))

## [1.4.0](https://github.com/xsitarcik/wrappers/compare/v1.3.0...v1.4.0) (2023-05-10)


### Features

* added picard rules for BedToIntervalList, MarkDuplicates and CreateSequenceDictionary ([dd384ea](https://github.com/xsitarcik/wrappers/commit/dd384eac4649dd88f10f57af745867bacd3a5ea2))

## [1.3.0](https://github.com/xsitarcik/wrappers/compare/v1.2.0...v1.3.0) (2023-05-10)


### Features

* added cutadapt for both single and paired ([0ce8037](https://github.com/xsitarcik/wrappers/commit/0ce803703acfaad3b45486b5e066019120ea8bdd))

## [1.2.0](https://github.com/xsitarcik/wrappers/compare/v1.1.0...v1.2.0) (2023-05-04)


### Features

* added basic samtools rules, faidx, index, sort and view ([a188213](https://github.com/xsitarcik/wrappers/commit/a1882135ea40ca48e9c55cdf605477a3a3543a27))
* added fastqc quality report ([878ca42](https://github.com/xsitarcik/wrappers/commit/878ca425aa9d2d32d09eb1df03815ff7c66c781e))
* added qualimap bamqc basic rule ([87b053c](https://github.com/xsitarcik/wrappers/commit/87b053ca1fee3e956edfcf9d7f14edea14e6911a))


### Bug Fixes

* move to new repo ([98925b2](https://github.com/xsitarcik/wrappers/commit/98925b23aa5042c0c26176623a15f99f0cf61ee1))

## [1.1.0](https://github.com/xsitarcik/wrappers/compare/v1.0.0...v1.1.0) (2023-05-03)


### Features

* added bwa mapping wrapper rules ([58e24d1](https://github.com/xsitarcik/wrappers/commit/58e24d18be0dd26e774829b11b3b4fb84c37f786))


### Bug Fixes

* bowtie2 wrapper to use correct wrapper repo instead of old one ([81bbc56](https://github.com/xsitarcik/wrappers/commit/81bbc56ed847756a1699e9be095c97b4fd82d993))
* bwa wrappers fixed to not use wildcards ([05ae11f](https://github.com/xsitarcik/wrappers/commit/05ae11f04f013aca0a14badc3d13e1ed07fced3b))
* bwa_mapping wrapper to use correct wrapper repo ([577142b](https://github.com/xsitarcik/wrappers/commit/577142b1fc345f0d996052a59408bcbf6dc22f78))

## 1.0.0 (2023-05-03)


### Features

* added bowtie2 mapping wrappers ([cde91e4](https://github.com/xsitarcik/wrappers/commit/cde91e403af67ae05f60b3ded0b61c0a6413853f))
