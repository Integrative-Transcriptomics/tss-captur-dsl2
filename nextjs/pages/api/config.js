// Define upload configuration
export const MAX_FILE_SIZE = (200) * 1024 * 1024; // 200MB
export const MAX_TOTAL_FILE_SIZE = (500) * 1024 * 1024; // 500MB
export const AVG_UPLOAD_SPEED = (1) / 8; // 0.125MB/s 
export const UPLOAD_BUFFER_TIME = 60 // in s
export const CHECK_INTERVAL = 5000; // in ms

// Handle configuration requests and responses
export default function configHandler(req, res) {
  try {
    // Check if the request method is GET
    if (req.method === 'GET') {
      // Respond with the configuration settings
      res.status(200).json({ maxFileSize: MAX_FILE_SIZE, maxTotalFileSize: MAX_TOTAL_FILE_SIZE, checkInterval: CHECK_INTERVAL });
    } else {
      // Respond with an error if the method is not allowed
      res.status(405).json({ error: 'Method not allowed' });
    }
  } catch (error) {
    // Handle any unexpected errors
    res.status(500).send('Internal Server Error');
  }
}

// Calculate the upload timeout based on the max total file size
export function calculateUploadTimeout() {
  // Calculate timeout in milliseconds
  return ((MAX_TOTAL_FILE_SIZE / (1024 * 1024) / AVG_UPLOAD_SPEED) + UPLOAD_BUFFER_TIME) * 1000; // in ms
}
