#!/bin/bash

#/*
# * @author Priyank Faldu <Priyank.Faldu@ed.ac.uk> <http://faldupriyank.com>
# *
# * Copyright 2019 The University of Edinburgh
# *
# * Licensed under the Apache License, Version 2.0 (the "License");
# * you may not use this file except in compliance with the License.
# * You may obtain a copy of the License at
# *
# *     http://www.apache.org/licenses/LICENSE-2.0
# *
# * Unless required by applicable law or agreed to in writing, software
# * distributed under the License is distributed on an "AS IS" BASIS,
# * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# * See the License for the specific language governing permissions and
# * limitations under the License.
# *
# */

DATASET="$1"

make -f ${DBG_ROOT}/graph-convert-utils/Makefile

#${DBG_ROOT}/graph-convert-utils/convert el bel ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert el vgr ${DATASET}
${DBG_ROOT}/graph-convert-utils/convert el cvgr ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert el csvgr ${DATASET}
${DBG_ROOT}/graph-convert-utils/convert el cintgr ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert bel vgr ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert bel cvgr ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert bel csvgr ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert bel cintgr ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert el print ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert bel print ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert vgr print ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert cvgr print ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert csvgr print ${DATASET}
#${DBG_ROOT}/graph-convert-utils/convert cintgr print ${DATASET}
