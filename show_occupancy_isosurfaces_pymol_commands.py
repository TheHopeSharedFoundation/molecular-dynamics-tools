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

obj_count = 0
obj_list =  cmd.get_object_list(selection='*3dBinned*')
for obj in obj_list:
	obj_count += 1
	cmd.create("%s_alpha_0.9to1.0"%obj_count,selection="(%s and resname ALP and b >0.9) and (%s and resname ALP and (b < 1.0 or b=1.0))"%(obj,obj))
	cmd.create("%s_alpha_0.8to0.9"%obj_count,selection="(%s and resname ALP and b >0.8) and (%s and resname ALP and (b < 0.9 or b=0.9))"%(obj,obj))
	cmd.create("%s_alpha_0.7to0.8"%obj_count,selection="(%s and resname ALP and b >0.7) and (%s and resname ALP and (b < 0.8 or b=0.8))"%(obj,obj))
	cmd.create("%s_alpha_0.6to0.7"%obj_count,selection="(%s and resname ALP and b >0.6) and (%s and resname ALP and (b < 0.7 or b=0.7))"%(obj,obj))
	cmd.create("%s_alpha_0.5to0.6"%obj_count,selection="(%s and resname ALP and b >0.5) and (%s and resname ALP and (b < 0.6 or b=0.6))"%(obj,obj))
	cmd.create("%s_alpha_0.4to0.5"%obj_count,selection="(%s and resname ALP and b >0.4) and (%s and resname ALP and (b < 0.5 or b=0.5))"%(obj,obj))
	cmd.create("%s_alpha_0.3to0.4"%obj_count,selection="(%s and resname ALP and b >0.3) and (%s and resname ALP and (b < 0.4 or b=0.4))"%(obj,obj))
	cmd.create("%s_alpha_0.2to0.3"%obj_count,selection="(%s and resname ALP and b >0.2) and (%s and resname ALP and (b < 0.3 or b=0.3))"%(obj,obj))
	cmd.create("%s_alpha_0.1to0.2"%obj_count,selection="(%s and resname ALP and b >0.1) and (%s and resname ALP and (b < 0.2 or b=0.2))"%(obj,obj))
	cmd.create("%s_alpha_0.0to0.1"%obj_count,selection="(%s and resname ALP and b >0.0) and (%s and resname ALP and (b < 0.1 or b=0.1))"%(obj,obj))
	cmd.group("%s_alpha"%obj_count,members="%s_alpha*"%obj_count)
	cmd.show(representation="mesh",selection="%s_alpha"%obj_count)
