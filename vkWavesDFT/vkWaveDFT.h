#pragma once

#include <vulkan/vulkan.h>
#include "vulkanexamplebase.h"
#include "VulkanDevice.hpp"
#include "VulkanBuffer.hpp"

#include "complex.h"
#include "vector.h"
#include "utils.h"

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
private:

	complex hTilde_0(int n_prime, int m_prime);

	float phillips(int n_prime, int m_prime);

	float dispersion(int n_prime, int m_prime);

	complex hTilde(float t, int n_prime, int m_prime);

	complex_vector_normal h_D_and_n(vector2 x, float t);

	void evaluateWaves(float t);

public:
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

	vkOceanWaveDFT();
	~vkOceanWaveDFT();

	//virtual void getEnabledFeatures();

	void buildCommandBuffers();

	void draw();

	void initMeshData();

	void setupVertexDescriptions();

	void setupDescriptorPool();

	void setupDescriptorSetLayout();

	void setupDescriptorSet();

	void preparePipelines();

	void prepareUniformBuffers();

	void updateUniformBuffers();

	void updateUniformBufferLight();

	void prepare();

	virtual void render();

	virtual void viewChanged();

};