#include "ImageSegment.h"


ImageSegment::ImageSegment(const Image& src, float sigma, float c, int min_size, int& num_ccs)
{
	Segment(src,sigma,c,min_size,num_ccs);
}


ImageSegment::~ImageSegment(void)
{
}


/*
sigma = 0.5, K = 500, min = 50
sigma = 0.5, K = 1000, min = 100.
*/
void ImageSegment::Segment(const Image& src, float sigma, float c, int min_size, int& num_ccs)
{
	int width = src.width;
	int height = src.height;

	Image smoothImg = ImageProcess::GausscianSeparateBlur(src,sigma);

	// build graph
	edge *edges = new edge[width*height*4];
	int num = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (x < width-1) {
				edges[num].a = y * width + x;
				edges[num].b = y * width + (x+1);
				edges[num].w = static_cast<float>(ImageExpandTools::EuclideanDistance(smoothImg.getPixel<float>(y,x),smoothImg.getPixel<float>(y,x+1)));
				num++;
			}

			if (y < height-1) {
				edges[num].a = y * width + x;
				edges[num].b = (y+1) * width + x;
				edges[num].w = static_cast<float>(ImageExpandTools::EuclideanDistance(smoothImg.getPixel<float>(y,x),smoothImg.getPixel<float>(y+1,x)));
				num++;
			}

			if ((x < width-1) && (y < height-1)) {
				edges[num].a = y * width + x;
				edges[num].b = (y+1) * width + (x+1);
				edges[num].w = static_cast<float>(ImageExpandTools::EuclideanDistance(smoothImg.getPixel<float>(y,x),smoothImg.getPixel<float>(y+1,x+1)));
				num++;
			}

			if ((x < width-1) && (y > 0)) {
				edges[num].a = y * width + x;
				edges[num].b = (y-1) * width + (x+1);
				edges[num].w = static_cast<float>(ImageExpandTools::EuclideanDistance(smoothImg.getPixel<float>(y,x),smoothImg.getPixel<float>(y-1,x+1)));
				num++;
			}
		}
	}


	// segment
	universe *uS = segmentGraph(width*height, num, edges, c);

	// post process small components
	for (int i = 0; i < num; i++) {
		int a = uS->find(edges[i].a);
		int b = uS->find(edges[i].b);
		if ((a != b) && ((uS->size(a) < min_size) || (uS->size(b) < min_size)))
			uS->join(a, b);
	}
	delete [] edges;
	num_ccs = uS->num_sets();


	Image output(width, height,src.imgType,3);

	// pick random colors for each component
	std::vector<Pixel<uchar>> colors(width*height);
	for (int i = 0; i < width*height; i++)
		colors[i] = colors[i].randColor();

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int comp = uS->find(y * width + x);
			ImageTools::setColorPixel(output,y, x, colors[comp]);
		}
	}  

	delete uS;

	string imgName = tools::fileNameFromTime(PATH,"ImageSeg",output.imgExt);
	output.save(imgName.c_str());
}



/*
 * Segment a graph
 *
 * Returns a disjoint-set forest representing the segmentation.
 *
 * num_vertices: number of vertices in graph.
 * num_edges: number of edges in graph
 * edges: array of edges.
 * c: constant for treshold function.
 */
universe* ImageSegment::segmentGraph(int num_vertices, int num_edges, edge *edges, float c) 
{ 
	// sort edges by weight
	std::sort(edges, edges + num_edges);

	// make a disjoint-set forest
	universe *u = new universe(num_vertices);

	// init thresholds
	float *threshold = new float[num_vertices];
	for (int i = 0; i < num_vertices; i++)
		threshold[i] = THRESHOLD(1,c);

	// for each edge, in non-decreasing weight order...
	for (int i = 0; i < num_edges; i++) {
		edge *pedge = &edges[i];

		// components conected by this edge
		int a = u->find(pedge->a);
		int b = u->find(pedge->b);
		if (a != b) {
			if ((pedge->w <= threshold[a]) && (pedge->w <= threshold[b])) {
					u->join(a, b);
					a = u->find(a);
					threshold[a] = pedge->w + THRESHOLD(u->size(a), c);
			}
		}
	}

	// free up
	delete threshold;
	return u;
}
