
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <gli/gli.hpp>

#include <vulkan/vulkan.h>
#include "vulkanexamplebase.h"
#include "VulkanDevice.hpp"
#include "VulkanBuffer.hpp"

#include "complex.h"
#include "vector.h"
#include "utils.h"

#define VERTEX_BUFFER_BIND_ID 0
#define ENABLE_VALIDATION false

#define M_PI       3.14159265358979323846   // pi

// Vertex layout for this example
struct vertexOcean {
	float pos[3];
	float normal[3];
	float htilde0[3];
	float htilde0Conj[3];
	float originalPos[3];
};

struct complex_vector_normal
{
	complex h; //wave height
	vector2 D; //displacement
	vector3 n; //normal
};


class vkOceanWaveDFT : public VulkanExampleBase
{
	complex hTilde_0(int n_prime, int m_prime) {
		complex r = gaussianRandomVariable();
		return r * sqrt(phillips(n_prime, m_prime) / 2.0f);
	}

	float phillips(int n_prime, int m_prime) {
		vector2 k(float(M_PI) * (2 * n_prime - N) / length,
			float(M_PI) * (2 * m_prime - N) / length);
		float k_length = k.length();
		if (k_length < 0.000001) return 0.0;

		float k_length2 = k_length * k_length;
		float k_length4 = k_length2 * k_length2;

		float k_dot_w = k.unit() * w.unit();
		float k_dot_w2 = k_dot_w * k_dot_w;

		float w_length = w.length();
		float L = w_length * w_length / g;
		float L2 = L * L;

		float damping = 0.001;
		float l2 = L2 * damping * damping;

		return A * exp(-1.0f / (k_length2 * L2)) / k_length4 * k_dot_w2 * exp(-k_length2 * l2);
	}

	float dispersion(int n_prime, int m_prime) {
		float w_0 = 2.0f * float(M_PI) / 200.0f;
		float kx = float(M_PI) * (2 * n_prime - N) / length;
		float kz = float(M_PI) * (2 * m_prime - N) / length;
		return floor(sqrt(g * sqrt(kx * kx + kz * kz)) / w_0) * w_0;
	}

	complex hTilde(float t, int n_prime, int m_prime) {
		int index = m_prime * Nplus1 + n_prime;

		complex htilde0(verticesOcean[index].htilde0[0], verticesOcean[index].htilde0[1]);
		complex htilde0mkconj(verticesOcean[index].htilde0Conj[0], verticesOcean[index].htilde0Conj[1]);

		float omegat = dispersion(n_prime, m_prime) * t;

		float cos_ = cos(omegat);
		float sin_ = sin(omegat);

		complex c0(cos_, sin_);
		complex c1(cos_, -sin_);

		complex res = htilde0 * c0 + htilde0mkconj * c1;

		return htilde0 * c0 + htilde0mkconj * c1;
	}

	complex_vector_normal h_D_and_n(vector2 x, float t) {
		complex h(0.0f, 0.0f);
		vector2 D(0.0f, 0.0f);
		vector3 n(0.0f, 0.0f, 0.0f);

		complex c, res, htilde_c;
		vector2 k;
		float kx, kz, k_length, k_dot_x;

		for (int m_prime = 0; m_prime < N; m_prime++) {
			kz = 2.0f * float(M_PI) * (m_prime - N / 2.0f) / length;
			for (int n_prime = 0; n_prime < N; n_prime++) {
				kx = 2.0f * float(M_PI) * (n_prime - N / 2.0f) / length;
				k = vector2(kx, kz);

				k_length = k.length();
				k_dot_x = k * x;

				c = complex(cos(k_dot_x), sin(k_dot_x));
				htilde_c = hTilde(t, n_prime, m_prime) * c;

				h = h + htilde_c;

				n = n + vector3(-kx * htilde_c.b, 0.0f, -kz * htilde_c.b);

				if (k_length < 0.000001) continue;
				D = D + vector2(kx / k_length * htilde_c.b, kz / k_length * htilde_c.b);
			}
		}

		n = (vector3(0.0f, 1.0f, 0.0f) - n).unit();

		complex_vector_normal cvn;
		cvn.h = h;
		cvn.D = D;
		cvn.n = n;
		return cvn;
	}

	void evaluateWaves(float t) {
		float lambda = -1.0;
		int index;
		vector2 x;
		vector2 d;
		complex_vector_normal h_d_and_n;
		for (int m_prime = 0; m_prime < N; m_prime++) 
		{
			for (int n_prime = 0; n_prime < N; n_prime++) 
			{
				index = m_prime * Nplus1 + n_prime;

				x = vector2(verticesOcean[index].pos[0], verticesOcean[index].pos[2]);

				h_d_and_n = h_D_and_n(x, t);

				verticesOcean[index].pos[1] = h_d_and_n.h.a;

				verticesOcean[index].pos[0] = verticesOcean[index].originalPos[0] + lambda * h_d_and_n.D.x;
				verticesOcean[index].pos[2] = verticesOcean[index].originalPos[1] + lambda * h_d_and_n.D.y;

				verticesOcean[index].normal[0] = h_d_and_n.n.x;
				verticesOcean[index].normal[1] = h_d_and_n.n.y;
				verticesOcean[index].normal[2] = h_d_and_n.n.z;

				if (n_prime == 0 && m_prime == 0) 
				{
					verticesOcean[index + N + Nplus1 * N].pos[1] = h_d_and_n.h.a;

					verticesOcean[index + N + Nplus1 * N].pos[0] = verticesOcean[index + N + Nplus1 * N].originalPos[0] + lambda * h_d_and_n.D.x;
					verticesOcean[index + N + Nplus1 * N].pos[2] = verticesOcean[index + N + Nplus1 * N].originalPos[2] + lambda * h_d_and_n.D.y;

					verticesOcean[index + N + Nplus1 * N].normal[0] = h_d_and_n.n.x;
					verticesOcean[index + N + Nplus1 * N].normal[1] = h_d_and_n.n.y;
					verticesOcean[index + N + Nplus1 * N].normal[2] = h_d_and_n.n.z;
				}
				if (n_prime == 0) 
				{
					verticesOcean[index + N].pos[1] = h_d_and_n.h.a;

					verticesOcean[index + N].pos[0] = verticesOcean[index + N].originalPos[0] + lambda * h_d_and_n.D.x;
					verticesOcean[index + N].pos[2] = verticesOcean[index + N].originalPos[2] + lambda * h_d_and_n.D.y;

					verticesOcean[index + N].normal[0] = h_d_and_n.n.x;
					verticesOcean[index + N].normal[1] = h_d_and_n.n.y;
					verticesOcean[index + N].normal[2] = h_d_and_n.n.z;
				}
				if (m_prime == 0) 
				{
					verticesOcean[index + Nplus1 * N].pos[1] = h_d_and_n.h.a;

					verticesOcean[index + Nplus1 * N].pos[0] = verticesOcean[index + Nplus1 * N].originalPos[0] + lambda * h_d_and_n.D.x;
					verticesOcean[index + Nplus1 * N].pos[2] = verticesOcean[index + Nplus1 * N].originalPos[2] + lambda * h_d_and_n.D.y;

					verticesOcean[index + Nplus1 * N].normal[0] = h_d_and_n.n.x;
					verticesOcean[index + Nplus1 * N].normal[1] = h_d_and_n.n.y;
					verticesOcean[index + Nplus1 * N].normal[2] = h_d_and_n.n.z;
				}
			}
		}

		size_t size = (Nplus1) * (Nplus1) * sizeof(vertexOcean);
		memcpy(oceanMappedMemory, verticesOcean, size);
	}
public:
	// Contains all Vulkan objects that are required to store and use a texture
	// Note that this repository contains a texture class (VulkanTexture.hpp) that encapsulates texture loading functionality in a class that is used in subsequent demos
	struct Texture {
		VkSampler sampler;
		VkImage image;
		VkImageLayout imageLayout;
		VkDeviceMemory deviceMemory;
		VkImageView view;
		uint32_t width, height;
		uint32_t mipLevels;
	} texture;

	struct {
		VkPipelineVertexInputStateCreateInfo inputState;
		std::vector<VkVertexInputBindingDescription> bindingDescriptions;
		std::vector<VkVertexInputAttributeDescription> attributeDescriptions;
	} vertices;

	vertexOcean* verticesOcean;
	unsigned int* indices;
	
	void* oceanMappedMemory;
	VkDeviceMemory oceanMemory;

	float g;				// gravity constant
	uint64_t N, Nplus1;				// dimension -- N should be a power of 2
	float A;				// phillips spectrum parameter -- affects heights of waves
	vector2 w;				// wind parameter
	float length;				// length parameter

	vks::Buffer vertexBuffer;
	vks::Buffer indexBuffer;
	uint32_t indexCount;

	vks::Buffer uniformBufferVS;

	struct {
		glm::mat4 projection;
		glm::mat4 model;
		glm::vec4 viewPos;
		glm::vec4 lightPos = glm::vec4(1000.0f, -100.0f, 1000.0f, 0.0f);
	} uboVS;

	struct {
		VkPipeline solid;
	} pipelines;

	VkPipelineLayout pipelineLayout;
	VkDescriptorSet descriptorSet;
	VkDescriptorSetLayout descriptorSetLayout;

	vkOceanWaveDFT() : VulkanExampleBase(ENABLE_VALIDATION)
	{
		zoom = -2.5f;
		rotation = { 0.0f, 15.0f, 0.0f };
		title = "ocean wave DFT";
		settings.overlay = true;

		N = 64;
		Nplus1 = N + 1;
		g = 9.81;
		A = 0.0005f;
		w = vector2(0.0f, 32.0f);
		length = 64;

		verticesOcean = new vertexOcean[Nplus1 * Nplus1];
		indices = new unsigned int[Nplus1 * Nplus1 * 10];

		oceanMappedMemory = nullptr;;
	}

	~vkOceanWaveDFT()
	{
		// Clean up used Vulkan resources 
		// Note : Inherited destructor cleans up resources stored in base class

		destroyTextureImage(texture);

		vkDestroyPipeline(device, pipelines.solid, nullptr);

		vkDestroyPipelineLayout(device, pipelineLayout, nullptr);
		vkDestroyDescriptorSetLayout(device, descriptorSetLayout, nullptr);

		//vkUnmapMemory(device, oceanMemory);
		//vkFreeMemory(device, oceanMemory, nullptr);

		if (verticesOcean)
		{
			delete[] verticesOcean;
		}

		if (indices)
		{
			delete[] indices;
		}

		vertexBuffer.destroy();
		indexBuffer.destroy();
		uniformBufferVS.destroy();
	}

	// Enable physical device features required for this example				
	virtual void getEnabledFeatures()
	{
		// Enable anisotropic filtering if supported
		if (deviceFeatures.samplerAnisotropy) {
			enabledFeatures.samplerAnisotropy = VK_TRUE;
		};
	}

	// Free all Vulkan resources used by a texture object
	void destroyTextureImage(Texture texture)
	{
		vkDestroyImageView(device, texture.view, nullptr);
		vkDestroyImage(device, texture.image, nullptr);
		vkDestroySampler(device, texture.sampler, nullptr);
		vkFreeMemory(device, texture.deviceMemory, nullptr);
	}

	void buildCommandBuffers()
	{
		VkCommandBufferBeginInfo cmdBufInfo = vks::initializers::commandBufferBeginInfo();

		VkClearValue clearValues[2];
		clearValues[0].color = defaultClearColor;
		clearValues[1].depthStencil = { 1.0f, 0 };

		VkRenderPassBeginInfo renderPassBeginInfo = vks::initializers::renderPassBeginInfo();
		renderPassBeginInfo.renderPass = renderPass;
		renderPassBeginInfo.renderArea.offset.x = 0;
		renderPassBeginInfo.renderArea.offset.y = 0;
		renderPassBeginInfo.renderArea.extent.width = width;
		renderPassBeginInfo.renderArea.extent.height = height;
		renderPassBeginInfo.clearValueCount = 2;
		renderPassBeginInfo.pClearValues = clearValues;

		for (int32_t i = 0; i < drawCmdBuffers.size(); ++i)
		{
			// Set target frame buffer
			renderPassBeginInfo.framebuffer = frameBuffers[i];

			VK_CHECK_RESULT(vkBeginCommandBuffer(drawCmdBuffers[i], &cmdBufInfo));

			vkCmdBeginRenderPass(drawCmdBuffers[i], &renderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

			VkViewport viewport = vks::initializers::viewport((float)width, (float)height, 0.0f, 1.0f);
			vkCmdSetViewport(drawCmdBuffers[i], 0, 1, &viewport);

			VkRect2D scissor = vks::initializers::rect2D(width, height, 0, 0);
			vkCmdSetScissor(drawCmdBuffers[i], 0, 1, &scissor);

			vkCmdBindDescriptorSets(drawCmdBuffers[i], VK_PIPELINE_BIND_POINT_GRAPHICS, pipelineLayout, 0, 1, &descriptorSet, 0, NULL);
			vkCmdBindPipeline(drawCmdBuffers[i], VK_PIPELINE_BIND_POINT_GRAPHICS, pipelines.solid);

			VkDeviceSize offsets[1] = { 0 };
			vkCmdBindVertexBuffers(drawCmdBuffers[i], VERTEX_BUFFER_BIND_ID, 1, &vertexBuffer.buffer, offsets);
			vkCmdBindIndexBuffer(drawCmdBuffers[i], indexBuffer.buffer, 0, VK_INDEX_TYPE_UINT32);

			vkCmdDrawIndexed(drawCmdBuffers[i], indexCount, 1, 0, 0, 0);

			drawUI(drawCmdBuffers[i]);

			vkCmdEndRenderPass(drawCmdBuffers[i]);

			VK_CHECK_RESULT(vkEndCommandBuffer(drawCmdBuffers[i]));
		}
	}

	void draw()
	{
		VulkanExampleBase::prepareFrame();

		// Command buffer to be sumitted to the queue
		submitInfo.commandBufferCount = 1;
		submitInfo.pCommandBuffers = &drawCmdBuffers[currentBuffer];

		// Submit to queue
		VK_CHECK_RESULT(vkQueueSubmit(queue, 1, &submitInfo, VK_NULL_HANDLE));

		VulkanExampleBase::submitFrame();
	}

	void initMeshData()
	{
		int index;

		complex htilde0, htilde0mk_conj;
		for (int m_prime = 0; m_prime < Nplus1; m_prime++) {
			for (int n_prime = 0; n_prime < Nplus1; n_prime++) {
				index = m_prime * Nplus1 + n_prime;

				htilde0 = hTilde_0(n_prime, m_prime);
				htilde0mk_conj = hTilde_0(-n_prime, -m_prime).conj();

				verticesOcean[index].htilde0[0] = htilde0.a;
				verticesOcean[index].htilde0[1] = htilde0.b;
				verticesOcean[index].htilde0Conj[0] = htilde0mk_conj.a;
				verticesOcean[index].htilde0Conj[1] = htilde0mk_conj.b;

				verticesOcean[index].originalPos[0] = verticesOcean[index].pos[0] = (n_prime - N / 2.0f) * length / N;
				verticesOcean[index].originalPos[1] = verticesOcean[index].pos[1] = 0.0f;
				verticesOcean[index].originalPos[2] = verticesOcean[index].pos[2] = (m_prime - N / 2.0f) * length / N;

				verticesOcean[index].normal[0] = 0.0f;
				verticesOcean[index].normal[1] = 1.0f;
				verticesOcean[index].normal[2] = 0.0f;
			}
		}

		indexCount = 0;
		for (int m_prime = 0; m_prime < N; m_prime++) 
		{
			for (int n_prime = 0; n_prime < N; n_prime++) 
			{
				index = m_prime * Nplus1 + n_prime;

				indices[indexCount++] = index;				// two triangles
				indices[indexCount++] = index + Nplus1;
				indices[indexCount++] = index + Nplus1 + 1;
				indices[indexCount++] = index;
				indices[indexCount++] = index + Nplus1 + 1;
				indices[indexCount++] = index + 1;
				
			}
		}

		size_t size = (Nplus1) * (Nplus1) * sizeof(vertexOcean);
		// Create buffers
		// For the sake of simplicity we won't stage the vertex data to the gpu memory
		// Vertex buffer
		VK_CHECK_RESULT(vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			size,
			&vertexBuffer.buffer,
			&oceanMemory,
			verticesOcean));
		// Index buffer
		VK_CHECK_RESULT(vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_INDEX_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			&indexBuffer,
			indexCount * sizeof(unsigned int),
			indices));

		VK_CHECK_RESULT(vkMapMemory(device, oceanMemory, 0, size, 0, &oceanMappedMemory));
	}

	void setupVertexDescriptions()
	{
		// Binding description
		vertices.bindingDescriptions.resize(1);
		vertices.bindingDescriptions[0] =
			vks::initializers::vertexInputBindingDescription(
				VERTEX_BUFFER_BIND_ID,
				sizeof(vertexOcean),
				VK_VERTEX_INPUT_RATE_VERTEX);

		// Attribute descriptions
		// Describes memory layout and shader positions
		vertices.attributeDescriptions.resize(2);
		// Location 0 : Position
		vertices.attributeDescriptions[0] =
			vks::initializers::vertexInputAttributeDescription(
				VERTEX_BUFFER_BIND_ID,
				0,
				VK_FORMAT_R32G32B32_SFLOAT,
				offsetof(vertexOcean, pos));
		// Location 1 : Texture coordinates
		vertices.attributeDescriptions[1] =
			vks::initializers::vertexInputAttributeDescription(
				VERTEX_BUFFER_BIND_ID,
				1,
				VK_FORMAT_R32G32B32_SFLOAT,
				offsetof(vertexOcean, normal));
		// Location 2 : Vertex normal
		/*vertices.attributeDescriptions[2] =
			vks::initializers::vertexInputAttributeDescription(
				VERTEX_BUFFER_BIND_ID,
				2,
				VK_FORMAT_R32G32B32_SFLOAT,
				offsetof(Vertex, normal));*/

		vertices.inputState = vks::initializers::pipelineVertexInputStateCreateInfo();
		vertices.inputState.vertexBindingDescriptionCount = static_cast<uint32_t>(vertices.bindingDescriptions.size());
		vertices.inputState.pVertexBindingDescriptions = vertices.bindingDescriptions.data();
		vertices.inputState.vertexAttributeDescriptionCount = static_cast<uint32_t>(vertices.attributeDescriptions.size());
		vertices.inputState.pVertexAttributeDescriptions = vertices.attributeDescriptions.data();
	}

	void setupDescriptorPool()
	{
		// Example uses one ubo and one image sampler
		std::vector<VkDescriptorPoolSize> poolSizes =
		{
			vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1)	
		};

		VkDescriptorPoolCreateInfo descriptorPoolInfo =
			vks::initializers::descriptorPoolCreateInfo(
				static_cast<uint32_t>(poolSizes.size()),
				poolSizes.data(),
				2);

		VK_CHECK_RESULT(vkCreateDescriptorPool(device, &descriptorPoolInfo, nullptr, &descriptorPool));
	}

	void setupDescriptorSetLayout()
	{
		std::vector<VkDescriptorSetLayoutBinding> setLayoutBindings =
		{
			// Binding 0 : Vertex shader uniform buffer
			vks::initializers::descriptorSetLayoutBinding(
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				VK_SHADER_STAGE_VERTEX_BIT,
				0)
		};

		VkDescriptorSetLayoutCreateInfo descriptorLayout = vks::initializers::descriptorSetLayoutCreateInfo(
				setLayoutBindings.data(),
				static_cast<uint32_t>(setLayoutBindings.size()));

		VK_CHECK_RESULT(vkCreateDescriptorSetLayout(device, &descriptorLayout, nullptr, &descriptorSetLayout));

		VkPipelineLayoutCreateInfo pPipelineLayoutCreateInfo = vks::initializers::pipelineLayoutCreateInfo(
				&descriptorSetLayout,
				1);

		VK_CHECK_RESULT(vkCreatePipelineLayout(device, &pPipelineLayoutCreateInfo, nullptr, &pipelineLayout));
	}

	void setupDescriptorSet()
	{
		VkDescriptorSetAllocateInfo allocInfo =
			vks::initializers::descriptorSetAllocateInfo(
				descriptorPool,
				&descriptorSetLayout,
				1);

		VK_CHECK_RESULT(vkAllocateDescriptorSets(device, &allocInfo, &descriptorSet));

		
		std::vector<VkWriteDescriptorSet> writeDescriptorSets =
		{
			// Binding 0 : Vertex shader uniform buffer
			vks::initializers::writeDescriptorSet(
				descriptorSet,
				VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
				0,
				&uniformBufferVS.descriptor),
			
		};

		vkUpdateDescriptorSets(device, static_cast<uint32_t>(writeDescriptorSets.size()), writeDescriptorSets.data(), 0, NULL);
	}

	void preparePipelines()
	{
		VkPipelineInputAssemblyStateCreateInfo inputAssemblyState =
			vks::initializers::pipelineInputAssemblyStateCreateInfo(
				VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST,
				0,
				VK_FALSE);

		VkPipelineRasterizationStateCreateInfo rasterizationState =
			vks::initializers::pipelineRasterizationStateCreateInfo(
				VK_POLYGON_MODE_FILL,
				VK_CULL_MODE_NONE,
				VK_FRONT_FACE_COUNTER_CLOCKWISE,
				0);

		VkPipelineColorBlendAttachmentState blendAttachmentState =
			vks::initializers::pipelineColorBlendAttachmentState(
				0xf,
				VK_FALSE);

		VkPipelineColorBlendStateCreateInfo colorBlendState =
			vks::initializers::pipelineColorBlendStateCreateInfo(
				1,
				&blendAttachmentState);

		VkPipelineDepthStencilStateCreateInfo depthStencilState =
			vks::initializers::pipelineDepthStencilStateCreateInfo(
				VK_TRUE,
				VK_TRUE,
				VK_COMPARE_OP_LESS_OR_EQUAL);

		VkPipelineViewportStateCreateInfo viewportState =
			vks::initializers::pipelineViewportStateCreateInfo(1, 1, 0);

		VkPipelineMultisampleStateCreateInfo multisampleState =
			vks::initializers::pipelineMultisampleStateCreateInfo(
				VK_SAMPLE_COUNT_1_BIT,
				0);

		std::vector<VkDynamicState> dynamicStateEnables = {
			VK_DYNAMIC_STATE_VIEWPORT,
			VK_DYNAMIC_STATE_SCISSOR
		};
		VkPipelineDynamicStateCreateInfo dynamicState =
			vks::initializers::pipelineDynamicStateCreateInfo(
				dynamicStateEnables.data(),
				static_cast<uint32_t>(dynamicStateEnables.size()),
				0);

		// Load shaders
		std::array<VkPipelineShaderStageCreateInfo, 2> shaderStages;

		shaderStages[0] = loadShader("shaders/oceanWaveDFT.vert.spv", VK_SHADER_STAGE_VERTEX_BIT);
		shaderStages[1] = loadShader("shaders/oceanWaveDFT.frag.spv", VK_SHADER_STAGE_FRAGMENT_BIT);

		VkGraphicsPipelineCreateInfo pipelineCreateInfo =
			vks::initializers::pipelineCreateInfo(
				pipelineLayout,
				renderPass,
				0);

		pipelineCreateInfo.pVertexInputState = &vertices.inputState;
		pipelineCreateInfo.pInputAssemblyState = &inputAssemblyState;
		pipelineCreateInfo.pRasterizationState = &rasterizationState;
		pipelineCreateInfo.pColorBlendState = &colorBlendState;
		pipelineCreateInfo.pMultisampleState = &multisampleState;
		pipelineCreateInfo.pViewportState = &viewportState;
		pipelineCreateInfo.pDepthStencilState = &depthStencilState;
		pipelineCreateInfo.pDynamicState = &dynamicState;
		pipelineCreateInfo.stageCount = static_cast<uint32_t>(shaderStages.size());
		pipelineCreateInfo.pStages = shaderStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(device, pipelineCache, 1, &pipelineCreateInfo, nullptr, &pipelines.solid));
	}

	// Prepare and initialize uniform buffer containing shader uniforms
	void prepareUniformBuffers()
	{
		// Vertex shader uniform buffer block
		VK_CHECK_RESULT(vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			&uniformBufferVS,
			sizeof(uboVS),
			&uboVS));

		updateUniformBuffers();
	}

	void updateUniformBuffers()
	{
		// Vertex shader
		uboVS.projection = glm::perspective(glm::radians(60.0f), (float)width / (float)height, 0.001f, 256.0f);
		glm::mat4 viewMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, zoom));

		uboVS.model = viewMatrix * glm::translate(glm::mat4(1.0f), cameraPos);
		uboVS.model = glm::rotate(uboVS.model, glm::radians(rotation.x), glm::vec3(1.0f, 0.0f, 0.0f));
		uboVS.model = glm::rotate(uboVS.model, glm::radians(rotation.y), glm::vec3(0.0f, 1.0f, 0.0f));
		uboVS.model = glm::rotate(uboVS.model, glm::radians(rotation.z), glm::vec3(0.0f, 0.0f, 1.0f));

		uboVS.viewPos = glm::vec4(0.0f, 0.0f, -zoom, 0.0f);

		VK_CHECK_RESULT(uniformBufferVS.map());
		memcpy(uniformBufferVS.mapped, &uboVS, sizeof(uboVS));
		//uniformBufferVS.unmap();
	}

	void updateUniformBufferLight()
	{
		// Environment
		uboVS.lightPos.x = sin(timer * 2.0f * float(M_PI)) * 1.5f;
		uboVS.lightPos.y = 0.0f;
		uboVS.lightPos.z = cos(timer * 2.0f * float(M_PI)) * 1.5f;
		memcpy(uniformBufferVS.mapped, &uboVS, sizeof(uboVS));
	}

	void prepare()
	{
		VulkanExampleBase::prepare();
		//loadTexture();
		initMeshData();
		setupVertexDescriptions();
		prepareUniformBuffers();
		setupDescriptorSetLayout();
		preparePipelines();
		setupDescriptorPool();
		setupDescriptorSet();
		buildCommandBuffers();
		prepared = true;
	}

	virtual void render()
	{
		if (!prepared)
			return;
		draw();
		evaluateWaves(frameTimer * 0.01f);
		updateUniformBufferLight();
	}

	virtual void viewChanged()
	{
		updateUniformBuffers();
	}

	virtual void OnUpdateUIOverlay(vks::UIOverlay* overlay)
	{
		/*if (overlay->header("Settings")) {
			if (overlay->sliderFloat("LOD bias", &uboVS.lodBias, 0.0f, (float)texture.mipLevels)) {
				updateUniformBuffers();
			}
		}*/
	}
};

VULKAN_EXAMPLE_MAIN()