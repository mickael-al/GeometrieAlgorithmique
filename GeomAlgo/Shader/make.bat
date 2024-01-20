"%VULKAN_SDK%\Bin\glslc.exe" shader.vert -o vert.spv
"%VULKAN_SDK%\Bin\glslc.exe" shader.frag -o frag.spv
"%VULKAN_SDK%\Bin\glslc.exe" shader_wireframe.frag -o shader_wireframe_fs.spv
"%VULKAN_SDK%\Bin\glslc.exe" water.vert -o water_vs.spv
"%VULKAN_SDK%\Bin\glslc.exe" water.frag -o water_fs.spv
"%VULKAN_SDK%\Bin\glslc.exe" shader_unlit.vert -o vert_unlit.spv
"%VULKAN_SDK%\Bin\glslc.exe" shader_unlit.frag -o frag_unlit.spv
"%VULKAN_SDK%\Bin\glslc.exe" skybox_vs.vert -o skybox_vs.spv
"%VULKAN_SDK%\Bin\glslc.exe" skybox_fs.frag -o skybox_fs.spv
"%VULKAN_SDK%\Bin\glslc.exe" Shadow.vert -o Shadow_vs.spv
"%VULKAN_SDK%\Bin\glslc.exe" Shadow.frag -o Shadow_fs.spv
"%VULKAN_SDK%\Bin\glslc.exe" shader_mirror.vert -o mirror_vs.spv
"%VULKAN_SDK%\Bin\glslc.exe" shader_mirror.frag -o mirror_fs.spv
"%VULKAN_SDK%\Bin\glslc.exe" particle.frag -o particle_fs.spv
"%VULKAN_SDK%\Bin\glslc.exe" particle.vert -o particle_vs.spv
"%VULKAN_SDK%\Bin\glslc.exe" compute_shader.comp -o compute_shader.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" bloom.comp -o bloom.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" tonemapping.comp -o tonemapping.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" bloom_mix.comp -o bloom_mix.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" chromatic_aberation.comp -o chromatic_aberation.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" depth_of_field.comp -o depth_of_field.comp.spv

"%VULKAN_SDK%\Bin\glslc.exe" nn\neural_network.comp -o nn\neural_network.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" nn\init_weight.comp -o nn\init_weight.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" nn\reset_score.comp -o nn\reset_score.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" nn\init_xor.comp -o nn\init_xor.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" nn\score_xor.comp -o nn\score_xor.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" nn\best_score.comp -o nn\best_score.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" nn\cmp_best.comp -o nn\cmp_best.comp.spv
"%VULKAN_SDK%\Bin\glslc.exe" nn\update_nn.comp -o nn\update_nn.comp.spv

"%VULKAN_SDK%\Bin\glslc.exe" real\character.frag -o real\character_fs.spv
"%VULKAN_SDK%\Bin\glslc.exe" real\character.vert -o real\character_vs.spv
"%VULKAN_SDK%\Bin\glslc.exe" real\pixel_water.frag -o real\pixel_water_fs.spv
"%VULKAN_SDK%\Bin\glslc.exe" real\pixel_water.vert -o real\pixel_water_vs.spv
"%VULKAN_SDK%\Bin\glslc.exe" real\noise.frag -o real\noise_fs.spv

pause

