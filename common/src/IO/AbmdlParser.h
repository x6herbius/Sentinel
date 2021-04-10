/*
 Copyright (C) 2010-2017 Kristian Duske

 This file is part of TrenchBroom.

 TrenchBroom is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 TrenchBroom is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with TrenchBroom. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "Assets/EntityModel_Forward.h"
#include "IO/EntityModelParser.h"

#include <vecmath/forward.h>
#include <vecmath/vec.h>
#include <vecmath/bbox.h>

#include <string>
#include <vector>

namespace TrenchBroom {
    namespace Assets {
        class Texture;
    }

    namespace IO {
        class BufferedReader;
        class Reader;
        class FileSystem;

        // This is not a nice class. We really need a proper MDL-handling framework to make it nice.
        class AbmdlParser : public EntityModelParser {
        private:
            using mat3x4f = vm::mat<float,3,4>;

            struct TextureInfo
            {
                std::string name;
                uint32_t flags = 0;
                int32_t width = 0;
                int32_t height = 0;
                int32_t index = 0;

                inline TextureInfo() = default;
                TextureInfo(BufferedReader& reader);
            };

            struct BodyPart
            {
                std::string name;
                size_t numModels = 0;
                int32_t base = 0;
                size_t modelIndex = 0;

                BodyPart(BufferedReader& reader);
            };

            struct SubModel
            {
                std::string name;
                int32_t unused = 0;
                float unused2;
                size_t numMeshes = 0;
                size_t meshIndex = 0;
                size_t numVertices = 0;
                size_t vertexInfoIndex = 0;
                size_t vertexIndex = 0;
                size_t numNormals = 0;
                size_t normalInfoIndex = 0;
                size_t normalIndex = 0;
                size_t blendVertexInfoIndex = 0;
                size_t blendNormalInfoIndex = 0;

                SubModel(BufferedReader& reader);
            };

            struct Mesh
            {
                size_t numTriangles = 0;
                size_t triangleIndex = 0;
                size_t skinRef = 0;
                size_t numNormals = 0;
                int32_t unused = 0;

                Mesh(BufferedReader& reader);
            };

            struct Bone
            {
                std::string name;
                int32_t parent = 0;
                int32_t unused = 0;
                int32_t controller[6] = { 0, 0, 0, 0, 0, 0 };
                float value[6] = { 0, 0, 0, 0, 0, 0 };
                float scale[6] = { 0, 0, 0, 0, 0, 0 };

                Bone(BufferedReader& reader);
            };

            struct TriangleVertex
            {
                size_t vertexIndex = 0;
                size_t normalIndex = 0;
                int32_t s = 0;
                int32_t t = 0;

                inline TriangleVertex() = default;
                TriangleVertex(BufferedReader& reader);
                void readFrom(BufferedReader& reader);
            };

            struct MdlIterator
            {
                size_t iteratorIndex;
                size_t bodyPartIndex;
                const BodyPart& bodyPart;
                const SubModel& firstSubModel;
                size_t meshIndex;
                const Mesh& mesh;

                std::string GenerateSurfaceName() const;
            };

            struct TriangleVertexTrio
            {
                const TriangleVertex* v[3];
            };

            struct VertexBatch
            {
                bool isTriangletrip = false;
                std::vector<Assets::EntityModelVertex> vertexList;
            };

            std::string m_name;
            const char* m_begin;
            const char* m_end;
            const FileSystem& m_fs;
        public:
            AbmdlParser(const std::string& name, const char* begin, const char* end, const FileSystem& fs);
        private:
            static vm::quatf anglesToQuaternion(float pitchRad, float yawRad, float rollRad);
            static mat3x4f quatAndOriginToMat(const vm::quatf& quat, const vm::vec3f& origin);
            static void concat3x4Matrices(const mat3x4f& a, const mat3x4f b, mat3x4f& out);
            static vm::vec3f transformVector(const vm::vec3f vec, const mat3x4f mat);

            std::unique_ptr<Assets::EntityModel> doInitializeModel(Logger& logger) override;
            void doLoadFrame(size_t frameIndex, Assets::EntityModel& model, Logger& logger) override;

            void loadTextures(Logger& logger, BufferedReader& reader, Assets::EntityModel& model, bool textureInfosOnly);
            void loadVerticesForSurfaces(BufferedReader& reader, Assets::EntityModel& model);
            void loadVerticesForSurface(BufferedReader& reader, const MdlIterator& it);
            void iterateMdlComponents(BufferedReader& reader, const std::function<void(const MdlIterator&)>& callback);

            void readBonesAndTransforms(BufferedReader& reader);
            void createBoneTransform(const Bone& bone);
            void readEmbeddedTexture(BufferedReader& reader, std::vector<Assets::Texture>* textureList, size_t textureIndex);
            void readTextureFromDisk(Logger& logger, BufferedReader& reader, std::vector<Assets::Texture>* textureList, size_t textureIndex);

            void cacheModelComponents(BufferedReader& reader, const SubModel& subModel);
            void getUnindexedVerticesFromMesh(BufferedReader& reader, std::vector<Assets::EntityModelVertex>& outVertices, const MdlIterator& it);
            void appendModelTriangle(std::vector<Assets::EntityModelVertex>& outVertices, const MdlIterator& it, const TriangleVertexTrio& verts);

            std::vector<Bone> m_bones;
            std::vector<mat3x4f> m_boneTransforms;
            std::vector<TextureInfo> m_textureInfos;

            size_t m_bodyPartIndexForCachedModelData = 0;
            std::vector<size_t> m_modelBoneIndices;
            std::vector<vm::vec3f> m_modelVertPositions;

            std::vector<VertexBatch> m_vertexBatches;
        };
    }
}
