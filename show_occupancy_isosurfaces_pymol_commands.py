############################################################################
# Copyright 2024 The Hope Shared Foundation                                #
#                                                                          #
# Licensed under the Apache License, Version 2.0 (the "License");          #
# you may not use this file except in compliance with the License.         #
# You may obtain a copy of the License at                                  #
#                                                                          #
# http://www.apache.org/licenses/LICENSE-2.0                               #
#                                                                          #
# Unless required by applicable law or agreed to in writing, software      #
# distributed under the License is distributed on an "AS IS" BASIS,        #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
# See the License for the specific language governing permissions and      #
# limitations under the License.                                           #
#                                                                          #
# The Hope Shared Foundation, September 2024                               #
# https://hopesharedfoundation.org/                                        #
# Towards Alleviating Suffering                                            #
############################################################################

cmd.create("alpha_0.9to1.0",selection="(resname ALP and b >0.9) and (resname ALP and (b < 1.0 or b=1.0))")
cmd.create("alpha_0.8to0.9",selection="(resname ALP and b >0.8) and (resname ALP and (b < 0.9 or b=0.9))")
cmd.create("alpha_0.7to0.8",selection="(resname ALP and b >0.7) and (resname ALP and (b < 0.8 or b=0.8))")
cmd.create("alpha_0.6to0.7",selection="(resname ALP and b >0.6) and (resname ALP and (b < 0.7 or b=0.7))")
cmd.create("alpha_0.5to0.6",selection="(resname ALP and b >0.5) and (resname ALP and (b < 0.6 or b=0.6))")
cmd.create("alpha_0.4to0.5",selection="(resname ALP and b >0.4) and (resname ALP and (b < 0.5 or b=0.5))")
cmd.create("alpha_0.3to0.4",selection="(resname ALP and b >0.3) and (resname ALP and (b < 0.4 or b=0.4))")
cmd.create("alpha_0.2to0.3",selection="(resname ALP and b >0.2) and (resname ALP and (b < 0.3 or b=0.3))")
cmd.create("alpha_0.1to0.2",selection="(resname ALP and b >0.1) and (resname ALP and (b < 0.2 or b=0.2))")
cmd.create("alpha_0.0to0.1",selection="(resname ALP and b >0.0) and (resname ALP and (b < 0.1 or b=0.1))")
cmd.group("alpha",members="alpha*")
cmd.show(representation="mesh",selection="alpha")
