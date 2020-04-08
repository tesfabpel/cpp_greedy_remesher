struct quad_t
{
    glm::ivec3 a, b, c, d;
};

static
size_t
to_index(const uint32_t u, 
    const uint32_t v, 
    const uint32_t sz_u)
{
    return v * sz_u
        + u;
}

static
std::vector<quad_t>
greedy_remesher(const glm::uvec3 &dims,
                const std::function<bool(glm::ivec3)> &map_fn)
{
	// https://github.com/AThilenius/greedy_meshing/blob/master/greedy_meshing.ts
	std::vector<quad_t> quads;

	// Sweep over 3-axes. `norm` is a number from [0, 3) which is used as an array
	// selection into a [x, y, z] vector3. Aka `norm` sweeps across the three
	// euclidean axises.
	for(size_t norm = 0u; norm < 3u; norm++)
	{
		// The tangent and bitangent axises. If you think of `norm` as the 'forward'
		// vector then `tan` is the 'right' vector, and `biTan` is the 'down'
		// vector. Again, each of these are just selectors into a [x, y, z] Vector3.
		const auto tan = (norm + 1u) % 3u;
		const auto biTan = (norm + 2u) % 3u;

		// The normal vector3, as a number triple (ideally you would replace this
		// with a THREE.Vector3 or something).
		glm::ivec3 normalVector(0);
		normalVector[norm] = 1;

		// Explained below.
		std::vector<bool> mask(dims[tan] * dims[biTan]);

		// Move through the volume in 2D 'slices' perpendicular to the
		// `normal_vector`. Iterate one more time for the "cap".
		for(size_t slice = 0u; slice <= dims[norm]; slice++)
		{
			// A 'voxel cursor' used to sample the 'slice' in the correct euclidean
			// plane.
			glm::uvec3 cursor(0);
			cursor[norm] = slice;

			// Compute the 2D mask of which voxels need to be tessellated.
			for(cursor[biTan] = 0; cursor[biTan] < dims[biTan]; ++cursor[biTan])
			{
				for(cursor[tan] = 0; cursor[tan] < dims[tan]; ++cursor[tan])
				{
					// The mask is set to true anywhere a voxel in the current 'slice`
					// differs from a voxel in the previous 'slice' in the
					// `-normalVector` direction. Aka anywhere a solid voxel face touches
					// a non-solid voxel. Note that this will cause sampling of
					// non-existent negative voxel coordinates, which
					// `getVoxelBoundsChecked` needs to handle by returning false.
					const glm::ivec3 curr(cursor);
					const glm::ivec3 vec(normalVector);

					const auto voxel_in_slice = map_fn(curr);
					const auto voxel_in_previous_slice = map_fn(curr - vec);

					const auto i = to_index(cursor[tan], cursor[biTan], dims[tan]);
					mask[i] = voxel_in_slice != voxel_in_previous_slice;
				}
			}

			// Generate mesh for mask using lexicographic ordering
			for(size_t y = 0; y < dims[biTan]; y++)
			{
				for(size_t x = 0; x < dims[tan];)
				{
					// If the mask isn't set, then just increment and continue (nothing to
					// tessellate).
					if(!mask[to_index(x, y, dims[tan])])
					{
						x++;
						continue;
					}

					// Compute the max-width of the combined quad going left-to-right.
					size_t width = 1;
					while(x + width < dims[tan]
						&& mask[to_index(x + width, y, dims[tan])])
					{
						width++;
					}

					// Compute max-height (extend the row `width` downward as much as
					// possible).
					size_t height = 1;
					for(; y + height < dims[biTan]; height++)
					{
						for(auto k = x; k < x + width; k++)
						{
							if(!mask[to_index(k, y + height, dims[tan])])
							{
								goto done_quad;
							}
						}
					}

				done_quad:
					// The base of the quad to add
					glm::ivec3 b(0);
					b[norm] = slice;
					b[tan] = x;
					b[biTan] = y;

					// The 'width' of the quad.
					glm::ivec3 du(0);
					du[tan] = width;

					// The 'height' of the quad.
					glm::ivec3 dv(0);
					dv[biTan] = height;

					quads.push_back({
						b,
						b + du,
						b + du + dv,
						b + dv,
					});

					// Clear the mask and increment x by the width of this quad.
					for(size_t l = 0; l < height; ++l)
					{
						for(size_t k = 0; k < width; ++k)
						{
							const auto i = to_index(x + k,
								y + l,
								dims[tan]);

							mask[i] = false;
						}
					}

					x += width;
				}
			}
		}
	}

	return quads;
}
