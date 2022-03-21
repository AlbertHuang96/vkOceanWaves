
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

#include "vkWaveDFT.h"


#define VERTEX_BUFFER_BIND_ID 0
#define ENABLE_VALIDATION false

#define M_PI       3.14159265358979323846   // pi



complex vkOceanWaveDFT::hTilde_0(int n_prime, int m_prime) 
{
	complex r = gaussianRandomVariable();
	return r * sqrt(phillips(n_prime, m_prime) / 2.0f);
}

float vkOceanWaveDFT::phillips(int n_prime, int m_prime) 
{
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

float vkOceanWaveDFT::dispersion(int n_prime, int m_prime) 
{
	float w_0 = 2.0f * float(M_PI) / 200.0f;
	float kx = float(M_PI) * (2 * n_prime - N) / length;
	float kz = float(M_PI) * (2 * m_prime - N) / length;
	return floor(sqrt(g * sqrt(kx * kx + kz * kz)) / w_0) * w_0;
}

complex vkOceanWaveDFT::hTilde(float t, int n_prime, int m_prime) {
	int index = m_prime * Nplus1 + n_prime;

	complex htilde0(oceanBuffer[index].htilde0.x, oceanBuffer[index].htilde0.y);
	complex htilde0mkconj(oceanBuffer[index].htilde0Conj.x, oceanBuffer[index].htilde0Conj.y);

	float omegat = dispersion(n_prime, m_prime) * t;

	float cos_ = cos(omegat);
	float sin_ = sin(omegat);

	complex c0(cos_, sin_);
	complex c1(cos_, -sin_);

	complex res = htilde0 * c0 + htilde0mkconj * c1;

	return htilde0 * c0 + htilde0mkconj * c1;
}

complex_vector_normal vkOceanWaveDFT::h_D_and_n(vector2 x, float t) {
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

void vkOceanWaveDFT::evaluateWaves(float t) {
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

			x = vector2(oceanBuffer[index].pos.x, oceanBuffer[index].pos.z);

			h_d_and_n = h_D_and_n(x, t);

			oceanBuffer[index].pos.y = h_d_and_n.h.a;

			oceanBuffer[index].pos.x = oceanBuffer[index].originalPos.x + lambda * h_d_and_n.D.x;
			oceanBuffer[index].pos.z = oceanBuffer[index].originalPos.z + lambda * h_d_and_n.D.y;

			oceanBuffer[index].normal.x = h_d_and_n.n.x;
			oceanBuffer[index].normal.y = h_d_and_n.n.y;
			oceanBuffer[index].normal.z = h_d_and_n.n.z;

			if (n_prime == 0 && m_prime == 0) 
			{
				oceanBuffer[index + N + Nplus1 * N].pos.y = h_d_and_n.h.a;

				oceanBuffer[index + N + Nplus1 * N].pos.x = oceanBuffer[index + N + Nplus1 * N].originalPos.x + lambda * h_d_and_n.D.x;
				oceanBuffer[index + N + Nplus1 * N].pos.z = oceanBuffer[index + N + Nplus1 * N].originalPos.z + lambda * h_d_and_n.D.y;

				oceanBuffer[index + N + Nplus1 * N].normal.x = h_d_and_n.n.x;
				oceanBuffer[index + N + Nplus1 * N].normal.y = h_d_and_n.n.y;
				oceanBuffer[index + N + Nplus1 * N].normal.z = h_d_and_n.n.z;
			}
			if (n_prime == 0) 
			{
				oceanBuffer[index + N].pos.y = h_d_and_n.h.a;

				oceanBuffer[index + N].pos.x = oceanBuffer[index + N].originalPos.x + lambda * h_d_and_n.D.x;
				oceanBuffer[index + N].pos.z = oceanBuffer[index + N].originalPos.z + lambda * h_d_and_n.D.y;

				oceanBuffer[index + N].normal.x = h_d_and_n.n.x;
				oceanBuffer[index + N].normal.y = h_d_and_n.n.y;
				oceanBuffer[index + N].normal.z = h_d_and_n.n.z;
			}
			if (m_prime == 0) 
			{
				oceanBuffer[index + Nplus1 * N].pos.y = h_d_and_n.h.a;

				oceanBuffer[index + Nplus1 * N].pos.x = oceanBuffer[index + Nplus1 * N].originalPos.x + lambda * h_d_and_n.D.x;
				oceanBuffer[index + Nplus1 * N].pos.z = oceanBuffer[index + Nplus1 * N].originalPos.z + lambda * h_d_and_n.D.y;

				oceanBuffer[index + Nplus1 * N].normal.x = h_d_and_n.n.x;
				oceanBuffer[index + Nplus1 * N].normal.y = h_d_and_n.n.y;
				oceanBuffer[index + Nplus1 * N].normal.z = h_d_and_n.n.z;
			}
		}
	}

	oceanWaves.size = (Nplus1) * (Nplus1) * sizeof(vertexOcean);
	memcpy(oceanWaves.mappedMemory, oceanBuffer.data(), oceanWaves.size);
}

void vkOceanWaveDFT::evaluateWavesFFT(float t) {
	float kx, kz, len, lambda = -1.0f;
	int index, index1;

	for (int m_prime = 0; m_prime < N; m_prime++) {
		kz = float(M_PI) * (2.0f * m_prime - N) / length;
		for (int n_prime = 0; n_prime < N; n_prime++) {
			kx = float(M_PI) * (2 * n_prime - N) / length;
			len = sqrt(kx * kx + kz * kz);
			index = m_prime * N + n_prime;

			h_tilde[index] = hTilde(t, n_prime, m_prime);
			h_tilde_slopex[index] = h_tilde[index] * complex(0, kx);
			h_tilde_slopez[index] = h_tilde[index] * complex(0, kz);
			if (len < 0.000001f) {
				h_tilde_dx[index] = complex(0.0f, 0.0f);
				h_tilde_dz[index] = complex(0.0f, 0.0f);
			}
			else {
				h_tilde_dx[index] = h_tilde[index] * complex(0, -kx / len);
				h_tilde_dz[index] = h_tilde[index] * complex(0, -kz / len);
			}
		}
	}

	for (int m_prime = 0; m_prime < N; m_prime++) {
		fft->fft(h_tilde, h_tilde, 1, m_prime * N);
		fft->fft(h_tilde_slopex, h_tilde_slopex, 1, m_prime * N);
		fft->fft(h_tilde_slopez, h_tilde_slopez, 1, m_prime * N);
		fft->fft(h_tilde_dx, h_tilde_dx, 1, m_prime * N);
		fft->fft(h_tilde_dz, h_tilde_dz, 1, m_prime * N);
	}
	for (int n_prime = 0; n_prime < N; n_prime++) {
		fft->fft(h_tilde, h_tilde, N, n_prime);
		fft->fft(h_tilde_slopex, h_tilde_slopex, N, n_prime);
		fft->fft(h_tilde_slopez, h_tilde_slopez, N, n_prime);
		fft->fft(h_tilde_dx, h_tilde_dx, N, n_prime);
		fft->fft(h_tilde_dz, h_tilde_dz, N, n_prime);
	}

	int sign;
	float signs[] = { 1.0f, -1.0f };
	vector3 n;
	for (int m_prime = 0; m_prime < N; m_prime++) {
		for (int n_prime = 0; n_prime < N; n_prime++) {
			index = m_prime * N + n_prime;		// index into h_tilde..
			index1 = m_prime * Nplus1 + n_prime;	// index into vertices

			sign = signs[(n_prime + m_prime) & 1];

			h_tilde[index] = h_tilde[index] * sign;

			// height
			oceanBuffer[index1].pos.y = h_tilde[index].a;

			// displacement
			h_tilde_dx[index] = h_tilde_dx[index] * sign;
			h_tilde_dz[index] = h_tilde_dz[index] * sign;
			oceanBuffer[index1].pos.x = oceanBuffer[index1].originalPos.x + h_tilde_dx[index].a * lambda;
			oceanBuffer[index1].pos.z = oceanBuffer[index1].originalPos.z + h_tilde_dz[index].a * lambda;

			// normal
			h_tilde_slopex[index] = h_tilde_slopex[index] * sign;
			h_tilde_slopez[index] = h_tilde_slopez[index] * sign;
			n = vector3(0.0f - h_tilde_slopex[index].a, 1.0f, 0.0f - h_tilde_slopez[index].a).unit();
			oceanBuffer[index1].normal.x = n.x;
			oceanBuffer[index1].normal.y = n.y;
			oceanBuffer[index1].normal.z = n.z;

			// for tiling
			if (n_prime == 0 && m_prime == 0) {
				oceanBuffer[index1 + N + Nplus1 * N].pos.y = h_tilde[index].a;

				oceanBuffer[index1 + N + Nplus1 * N].pos.x = oceanBuffer[index1 + N + Nplus1 * N].originalPos.x + h_tilde_dx[index].a * lambda;
				oceanBuffer[index1 + N + Nplus1 * N].pos.z = oceanBuffer[index1 + N + Nplus1 * N].originalPos.z + h_tilde_dz[index].a * lambda;

				oceanBuffer[index1 + N + Nplus1 * N].normal.x = n.x;
				oceanBuffer[index1 + N + Nplus1 * N].normal.y = n.y;
				oceanBuffer[index1 + N + Nplus1 * N].normal.z = n.z;
			}
			if (n_prime == 0) {
				oceanBuffer[index1 + N].pos.y = h_tilde[index].a;

				oceanBuffer[index1 + N].pos.x = oceanBuffer[index1 + N].originalPos.x + h_tilde_dx[index].a * lambda;
				oceanBuffer[index1 + N].pos.z = oceanBuffer[index1 + N].originalPos.z + h_tilde_dz[index].a * lambda;

				oceanBuffer[index1 + N].normal.x = n.x;
				oceanBuffer[index1 + N].normal.y = n.y;
				oceanBuffer[index1 + N].normal.z = n.z;
			}
			if (m_prime == 0) {
				oceanBuffer[index1 + Nplus1 * N].pos.y = h_tilde[index].a;

				oceanBuffer[index1 + Nplus1 * N].pos.x = oceanBuffer[index1 + Nplus1 * N].originalPos.x + h_tilde_dx[index].a * lambda;
				oceanBuffer[index1 + Nplus1 * N].pos.z = oceanBuffer[index1 + Nplus1 * N].originalPos.z + h_tilde_dz[index].a * lambda;

				oceanBuffer[index1 + Nplus1 * N].normal.x = n.x;
				oceanBuffer[index1 + Nplus1 * N].normal.y = n.y;
				oceanBuffer[index1 + Nplus1 * N].normal.z = n.z;
			}
		}
	}

	oceanWaves.size = (Nplus1) * (Nplus1) * sizeof(vertexOcean);
	memcpy(oceanWaves.mappedMemory, oceanBuffer.data(), oceanWaves.size);
}



vkOceanWaveDFT::vkOceanWaveDFT() : VulkanExampleBase(ENABLE_VALIDATION)
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

	h_tilde = new complex[N * N];
	h_tilde_slopex = new complex[N * N];
	h_tilde_slopez = new complex[N * N];
	h_tilde_dx = new complex[N * N];
	h_tilde_dz = new complex[N * N];
	fft = new cFFT(N);

	//verticesOcean = new vertexOcean[Nplus1 * Nplus1];
	//indices = new unsigned int[Nplus1 * Nplus1 * 10];

	//oceanMappedMemory = nullptr;;
}

vkOceanWaveDFT::~vkOceanWaveDFT()
{
	// Clean up used Vulkan resources 
	// Note : Inherited destructor cleans up resources stored in base class


	vkDestroyPipeline(device, pipelines.solid, nullptr);

	vkDestroyPipelineLayout(device, pipelineLayout, nullptr);
	vkDestroyDescriptorSetLayout(device, descriptorSetLayout, nullptr);

	if (h_tilde)		delete[] h_tilde;
	if (h_tilde_slopex)	delete[] h_tilde_slopex;
	if (h_tilde_slopez)	delete[] h_tilde_slopez;
	if (h_tilde_dx)		delete[] h_tilde_dx;
	if (h_tilde_dz)		delete[] h_tilde_dz;
	if (fft)		delete fft;

	/*if (verticesOcean)
	{
		delete[] verticesOcean;
	}*/

	/*if (indices)
	{
		delete[] indices;
	}*/

	vkUnmapMemory(device, oceanWaves.memory);
	vkDestroyBuffer(device, oceanWaves.buffer, nullptr);
	vkFreeMemory(device, oceanWaves.memory, nullptr);

	indexBuffer.destroy();
	uniformBufferVS.destroy();
}



void vkOceanWaveDFT::buildCommandBuffers()
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
		vkCmdBindVertexBuffers(drawCmdBuffers[i], VERTEX_BUFFER_BIND_ID, 1, &oceanWaves.buffer, offsets);
		vkCmdBindIndexBuffer(drawCmdBuffers[i], indexBuffer.buffer, 0, VK_INDEX_TYPE_UINT32);

		vkCmdDrawIndexed(drawCmdBuffers[i], indexCount, 1, 0, 0, 0);
		//vkCmdDraw(drawCmdBuffers[i], (Nplus1) * (Nplus1), 1, 0, 0);

		drawUI(drawCmdBuffers[i]);

		vkCmdEndRenderPass(drawCmdBuffers[i]);

		VK_CHECK_RESULT(vkEndCommandBuffer(drawCmdBuffers[i]));
	}
}

void vkOceanWaveDFT::draw()
{
	VulkanExampleBase::prepareFrame();

	// Command buffer to be sumitted to the queue
	submitInfo.commandBufferCount = 1;
	submitInfo.pCommandBuffers = &drawCmdBuffers[currentBuffer];

	// Submit to queue
	VK_CHECK_RESULT(vkQueueSubmit(queue, 1, &submitInfo, VK_NULL_HANDLE));

	VulkanExampleBase::submitFrame();
}

void vkOceanWaveDFT::initMeshData()
{
	oceanBuffer.resize(Nplus1 * Nplus1);

	int index;
	
	complex htilde0, htilde0mk_conj;
	for (int m_prime = 0; m_prime < Nplus1; m_prime++) {
		for (int n_prime = 0; n_prime < Nplus1; n_prime++) {
			index = m_prime * Nplus1 + n_prime;

			htilde0 = hTilde_0(n_prime, m_prime);
			htilde0mk_conj = hTilde_0(-n_prime, -m_prime).conj();

			oceanBuffer[index].htilde0.x = htilde0.a;
			oceanBuffer[index].htilde0.y = htilde0.b;
			oceanBuffer[index].htilde0Conj.x = htilde0mk_conj.a;
			oceanBuffer[index].htilde0Conj.y = htilde0mk_conj.b;

			oceanBuffer[index].originalPos.x = oceanBuffer[index].pos.x = (n_prime - N / 2.0f) * length / N;
			oceanBuffer[index].originalPos.y = oceanBuffer[index].pos.y = 0.0f;
			oceanBuffer[index].originalPos.z = oceanBuffer[index].pos.z = (m_prime - N / 2.0f) * length / N;

			oceanBuffer[index].normal.x = 0.0f;
			oceanBuffer[index].normal.y = 1.0f;
			oceanBuffer[index].normal.z = 0.0f;
		}
	}

	indices.resize(Nplus1 * Nplus1 * 10);
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

	oceanWaves.size = (Nplus1) * (Nplus1) * sizeof(vertexOcean);
	// Create buffers
	// For the sake of simplicity we won't stage the vertex data to the gpu memory
	// Vertex buffer
	VK_CHECK_RESULT(vulkanDevice->createBuffer(
		VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
		VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
		oceanWaves.size,
		&oceanWaves.buffer,
		&oceanWaves.memory,
		oceanBuffer.data()));
	// Index buffer
	VK_CHECK_RESULT(vulkanDevice->createBuffer(
		VK_BUFFER_USAGE_INDEX_BUFFER_BIT,
		VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
		&indexBuffer,
		indexCount * sizeof(uint32_t),
		indices.data()));

	VK_CHECK_RESULT(vkMapMemory(device, oceanWaves.memory, 0, oceanWaves.size, 0, &oceanWaves.mappedMemory));
}

void vkOceanWaveDFT::setupVertexDescriptions()
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
	vertices.attributeDescriptions[0] = vks::initializers::vertexInputAttributeDescription(VERTEX_BUFFER_BIND_ID, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(vertexOcean, pos));
	// Location 1 : Texture coordinates
	vertices.attributeDescriptions[1] = vks::initializers::vertexInputAttributeDescription(VERTEX_BUFFER_BIND_ID, 1, VK_FORMAT_R32G32B32_SFLOAT, offsetof(vertexOcean, normal));

	vertices.inputState = vks::initializers::pipelineVertexInputStateCreateInfo();
	vertices.inputState.vertexBindingDescriptionCount = static_cast<uint32_t>(vertices.bindingDescriptions.size());
	vertices.inputState.pVertexBindingDescriptions = vertices.bindingDescriptions.data();
	vertices.inputState.vertexAttributeDescriptionCount = static_cast<uint32_t>(vertices.attributeDescriptions.size());
	vertices.inputState.pVertexAttributeDescriptions = vertices.attributeDescriptions.data();
}

void vkOceanWaveDFT::setupDescriptorPool()
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

void vkOceanWaveDFT::setupDescriptorSetLayout()
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

void vkOceanWaveDFT::setupDescriptorSet()
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

void vkOceanWaveDFT::preparePipelines()
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
void vkOceanWaveDFT::prepareUniformBuffers()
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

void vkOceanWaveDFT::updateUniformBuffers()
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

void vkOceanWaveDFT::updateUniformBufferLight()
{
	// Environment
	uboVS.lightPos.x = sin(timer * 2.0f * float(M_PI)) * 1.5f;
	uboVS.lightPos.y = 0.0f;
	uboVS.lightPos.z = cos(timer * 2.0f * float(M_PI)) * 1.5f;
	memcpy(uniformBufferVS.mapped, &uboVS, sizeof(uboVS));
}

void vkOceanWaveDFT::prepare()
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

void vkOceanWaveDFT::render()
{
	if (!prepared)
		return;
	draw();
	//evaluateWaves(timer); // the frameTimer factor can be adjusted by multiplied by a scalar
	evaluateWavesFFT(timer * 3.0f);
	updateUniformBufferLight();
}

void vkOceanWaveDFT::viewChanged()
{
	updateUniformBuffers();
}


VULKAN_EXAMPLE_MAIN()